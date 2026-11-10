//! report_run — Markdown report for one individual from a simple_ga run.
//!
//! Usage:
//!   report_run --id <individual_id>  [--data-dir <path>]  [--seed <u64>]  [--out <file>]
//!
//! Reads  <data-dir>/phylogeny.json          (config, target, genome, scores)
//!   and  <data-dir>/reaction_logs.jsonl     (original reaction sequence)
//! Re-runs the simulation from scratch, computes all similarity metrics,
//! and writes a self-contained Markdown file to
//!   <data-dir>/markdown/individual_<id>.md  (or --out path).

use clap::Parser;
use rand::{SeedableRng, rngs::StdRng};
use serde_json;
use ssa_model::Gene8::Gene8;
use ssa_model::molecules::Molecule;
use ssa_model::reactions::Reaction;
use ssa_model::reactor::Reactor;
use std::collections::HashMap;
use std::fmt::Write;
const FORM: u8            = 1;
pub fn build_target() -> Molecule {
    let mut mol = Molecule::new();
    mol.name = "target_dimer".to_string();
    mol.atoms.insert(0, (FORM, 1));
    mol.atoms.insert(1, (FORM, 2));
    mol.atoms.insert(2, (FORM, 2));
    mol.atoms.insert(3, (FORM, 2));
    mol.atoms.insert(4, (FORM, 1));
    mol.links.insert(0, (0, 1));
    mol.links.insert(1, (1, 2));
    mol.links.insert(2, (2, 3));
    mol.links.insert(3, (3, 4));
    mol
}
// ── CLI ───────────────────────────────────────────────────────────────────────
#[derive(Parser)]
#[command(name = "report_run",
          about = "Markdown report for one individual from a simple_ga run")]
struct Cli {
    /// Individual ID to report
    #[arg(long)]
    id: usize,

    /// Directory produced by simple_ga (contains phylogeny.json etc.)
    #[arg(long, default_value = "results/simple_ga")]
    data_dir: String,

    /// RNG seed for re-run simulation (0 = use the seed stored in run config)
    #[arg(long, default_value = "0")]
    seed: u64,

    /// Output path (default: <data-dir>/markdown/individual_<id>.md)
    #[arg(long)]
    out: Option<String>,
}

// ── Reaction type label ───────────────────────────────────────────────────────
fn rxn_type(r: &Reaction) -> &'static str {
    if r.is_production()  { return "production"; }
    if r.is_degradation() { return "degradation"; }
    match (r.is_substrate_contact(), r.is_product_contact()) {
        (true,  false) => "contact→bond",
        (true,  true)  => "contact",
        (false, true)  => "link→split",
        (false, false) => "link→state",
    }
}

// ── Snapshot formatting helpers ───────────────────────────────────────────────
/// atom key (form<<8 | state) → "a0", "b3", …
fn key_to_sym(k: u16) -> String {
    let form  = ((k >> 8) & 0xff) as u8;
    let state = (k & 0xff) as u8;
    format!("{}{}", (form + 96) as char, state)
}

/// u32 link key → normalised "aX-aY" (lex smaller first)
fn normalize_link_key(k: u32) -> String {
    let a = key_to_sym((k >> 16) as u16);
    let b = key_to_sym((k & 0xffff) as u16);
    if a <= b { format!("{a}-{b}") } else { format!("{b}-{a}") }
}

/// Format a reactor snapshot as compact (atoms_str, links_str).
/// atoms_str: "a0:8 a1:1 a2:1"   links_str: "a1-a2:2" or "—"
fn format_snapshot(
    atom_snap: &HashMap<u16, usize>,
    link_snap: &HashMap<u32, usize>,
) -> (String, String) {
    let mut atom_parts: Vec<(String, usize)> = atom_snap.iter()
        .map(|(&k, &v)| (key_to_sym(k), v))
        .collect();
    atom_parts.sort_by_key(|(s, _)| s.clone());
    let atoms_str = atom_parts.iter()
        .map(|(s, v)| format!("{s}:{v}"))
        .collect::<Vec<_>>().join(" ");

    let mut link_map: HashMap<String, usize> = HashMap::new();
    for (&k, &v) in link_snap.iter() {
        *link_map.entry(normalize_link_key(k)).or_insert(0) += v;
    }
    let mut link_parts: Vec<_> = link_map.iter().collect();
    link_parts.sort_by_key(|(k, _)| k.clone());
    let links_str = if link_parts.is_empty() {
        "—".to_string()
    } else {
        link_parts.iter().map(|(s, v)| format!("{s}:{v}")).collect::<Vec<_>>().join(" ")
    };
    (atoms_str, links_str)
}

// ── Target reconstruction from phylogeny.json ─────────────────────────────────
fn target_from_json(v: &serde_json::Value) -> Molecule {
    let mut mol = Molecule::new();
    mol.name = v["name"].as_str().unwrap_or("target").to_string();
    if let Some(atoms) = v["atoms"].as_array() {
        for a in atoms {
            mol.atoms.insert(
                a["id"].as_u64().unwrap_or(0) as u32,
                (a["form"].as_u64().unwrap_or(0) as u8,
                 a["state"].as_u64().unwrap_or(0) as u8),
            );
        }
    }
    if let Some(links) = v["links"].as_array() {
        for (i, l) in links.iter().enumerate() {
            mol.links.insert(i as u32,
                (l["a"].as_u64().unwrap_or(0) as u32,
                 l["b"].as_u64().unwrap_or(0) as u32));
        }
    }
    mol
}

// ── Score / similarity ────────────────────────────────────────────────────────
struct ScoreDetail {
    n_linked:    usize,
    best_comp:   f64,
    best_size:   f64,
    iso_found:   bool,
    score_linked: f64,
    score_comp:   f64,
    score_size:   f64,
    score_iso:    f64,
    score_total:  f64,
}

fn score_full(reactor: &Reactor, target: &Molecule, iso_max: usize) -> ScoreDetail {
    let mols   = reactor.find_molecules();
    let linked: Vec<&Molecule> = mols.iter().filter(|m| m.size() > 1).collect();
    let n_linked = linked.len();
    let t_size   = target.size() as f64;

    let best_comp = linked.iter().map(|m| m.composition_similarity(target))
        .fold(0.0_f64, f64::max);
    let best_size = if linked.is_empty() { 0.0 } else {
        linked.iter()
            .map(|m| 1.0 / (1.0 + (m.size() as f64 - t_size).abs() / t_size.max(1.0)))
            .fold(0.0_f64, f64::max)
    };
    let iso_found = if target.size() <= iso_max {
        linked.iter().filter(|m| m.size() == target.size())
            .any(|m| m.composition_similarity(target) > 0.8 && m.is_isomorphic_to(target))
    } else { false };

    let sl = n_linked as f64 / 2.0;
    let sc = best_comp * 10.0;
    let ss = best_size * 5.0;
    let si = if iso_found { 100.0 } else { 0.0 };
    ScoreDetail { n_linked, best_comp, best_size, iso_found,
        score_linked: sl, score_comp: sc, score_size: ss, score_iso: si,
        score_total: sl + sc + ss + si }
}

// ── Markdown builder ──────────────────────────────────────────────────────────
fn generate_markdown(
    ind_id:      usize,
    node:        &serde_json::Value,
    gene:        &Gene8,
    target:      &Molecule,
    initial_key: u16,
    n_atoms:     usize,
    tmax:        f64,
    iso_max:     usize,
    rerun_seed:  u64,
    orig_log:    Option<&serde_json::Value>,
) -> String {
    let mut md = String::new();

    // ── Title ─────────────────────────────────────────────────────────────────
    let _ = writeln!(md, "# Individual Report — id {}", ind_id);
    let _ = writeln!(md, "> Generated by `report_run`  \n");

    // ── Metadata ─────────────────────────────────────────────────────────────
    let _ = writeln!(md, "## Individual Info\n");
    let _ = writeln!(md, "| Field | Value |");
    let _ = writeln!(md, "|-------|-------|");
    let _ = writeln!(md, "| ID | {} |", ind_id);
    let _ = writeln!(md, "| Generation | {} |", node["generation"]);
    let _ = writeln!(md, "| Position in pop | {} |", node["position"]);
    let _ = writeln!(md, "| Score (original) | **{:.4}** |",
        node["score"].as_f64().unwrap_or(0.0));
    let _ = writeln!(md, "| Decoded reactions | {} |",
        node["n_decoded_reactions"]);
    let _ = writeln!(md, "| Elite | {} |",
        node["is_elite"].as_bool().unwrap_or(false));
    let p1 = if node["parent1_id"].is_null() { "—".to_string() }
             else { node["parent1_id"].to_string() };
    let p2 = if node["parent2_id"].is_null() { "—".to_string() }
             else { node["parent2_id"].to_string() };
    let _ = writeln!(md, "| Parent 1 | {} |", p1);
    let _ = writeln!(md, "| Parent 2 | {} |", p2);
    let _ = writeln!(md);

    // ── Genome ────────────────────────────────────────────────────────────────
    let genome_str = gene.to_string();
    let _ = writeln!(md, "## Genome\n");
    let _ = writeln!(md, "- **Length**: {} tokens", genome_str.len());

    // Count E/F patterns
    let e_count = genome_str.chars().filter(|&c| c == 'E').count();
    let f_count = genome_str.chars().filter(|&c| c == 'F').count();
    let _ = writeln!(md, "- **START tokens (E)**: {}  **END tokens (F)**: {}",
        e_count, f_count);
    let _ = writeln!(md);
    let _ = writeln!(md, "```\n{}\n```\n", genome_str);

    // Annotated with markers highlighted
    let _ = writeln!(md, "<details><summary>Annotated (E=pattern start, F=end)</summary>\n");
    let annotated: String = genome_str.chars()
        .map(|c| match c { 'E' => "[E]".to_string(), 'F' => "[F]".to_string(),
                           other => other.to_string() })
        .collect();
    let _ = writeln!(md, "```\n{}\n```\n", annotated);
    let _ = writeln!(md, "</details>\n");

    // ── Decoded reactions ─────────────────────────────────────────────────────
    let reactions = gene.to_reactions();
    let _ = writeln!(md, "## Decoded Reactions ({} valid)\n", reactions.len());

    if reactions.is_empty() {
        let _ = writeln!(md, "*No valid reactions decoded from this genome.*\n");
    } else {
        let _ = writeln!(md, "| # | Reaction string | Type | Speed |");
        let _ = writeln!(md, "|---|-----------------|------|-------|");
        for (i, r) in reactions.iter().enumerate() {
            let _ = writeln!(md, "| {} | `{}` | {} | {:.3} |",
                i, r.to_string(), rxn_type(r), r.reaction_speed());
        }
        let _ = writeln!(md);
    }

    // ── Original run (from reaction_logs.jsonl) ───────────────────────────────
    if let Some(log) = orig_log {
        let total = log["total_reactions"].as_u64().unwrap_or(0);
        let _ = writeln!(md, "## Original Run (from log)\n");
        let _ = writeln!(md, "**Total reactions fired**: {}\n", total);

        if let Some(seq) = log["reaction_sequence"].as_array() {
            if seq.is_empty() {
                let _ = writeln!(md, "*No reactions fired in original run.*\n");
            } else {
                let show = seq.len().min(100);
                if seq.len() > 100 {
                    let _ = writeln!(md,
                        "*Showing first 100 of {} events — see `reaction_logs.jsonl` for full sequence.*\n",
                        seq.len());
                } else {
                    let _ = writeln!(md, "### Reaction Sequence\n");
                }
                let _ = writeln!(md, "| Step | Time | Reaction | Atoms (after) | Bonds (after) |");
                let _ = writeln!(md, "|------|------|----------|---------------|---------------|");
                for (i, ev) in seq.iter().take(show).enumerate() {
                    // Parse atoms/links maps if present (new format); fall back to "—"
                    let atoms_str = if let Some(obj) = ev["atoms"].as_object() {
                        let mut parts: Vec<_> = obj.iter()
                            .map(|(k, v)| (k.clone(), v.as_u64().unwrap_or(0) as usize))
                            .collect();
                        parts.sort_by_key(|(k, _)| k.clone());
                        parts.iter().map(|(k, v)| format!("{k}:{v}")).collect::<Vec<_>>().join(" ")
                    } else { "—".to_string() };
                    let links_str = if let Some(obj) = ev["links"].as_object() {
                        if obj.is_empty() { "—".to_string() } else {
                            let mut parts: Vec<_> = obj.iter()
                                .map(|(k, v)| (k.clone(), v.as_u64().unwrap_or(0) as usize))
                                .collect();
                            parts.sort_by_key(|(k, _)| k.clone());
                            parts.iter().map(|(k, v)| format!("{k}:{v}")).collect::<Vec<_>>().join(" ")
                        }
                    } else { "—".to_string() };
                    let _ = writeln!(md, "| {} | {:.6} | `{}` | {} | {} |",
                        i + 1,
                        ev["time"].as_f64().unwrap_or(0.0),
                        ev["reaction_str"].as_str().unwrap_or("?"),
                        atoms_str, links_str);
                }
                let _ = writeln!(md);
            }
        }

        if let Some(fmols) = log["final_molecules"].as_array() {
            let _ = writeln!(md, "### Final Molecules (original)\n");
            let _ = writeln!(md, "| Size | Atoms (form,state) | Bonds |");
            let _ = writeln!(md, "|------|--------------------|-------|");
            let mut shown = std::collections::HashMap::<String, usize>::new();
            for m in fmols {
                let atoms_arr = m["atoms"].as_array();
                let atom_str = atoms_arr
                    .map(|a| a.iter()
                        .map(|v| {
                            let r = vec![];
                            let pair = v.as_array().unwrap_or(&r);
                            format!("({},{})",
                                pair.get(0).and_then(|x| x.as_u64()).unwrap_or(0),
                                pair.get(1).and_then(|x| x.as_u64()).unwrap_or(0))
                        })
                        .collect::<Vec<_>>().join(" ")
                    )
                    .unwrap_or_default();
                let key = format!("{}-{}", m["size"], atom_str);
                *shown.entry(key.clone()).or_insert(0) += 1;
            }
            let mut rows: Vec<(String, usize)> = shown.into_iter().collect();
            rows.sort_by_key(|(k, _)| k.clone());
            for (desc, count) in &rows {
                let parts: Vec<&str> = desc.splitn(2, '-').collect();
                let sz = parts.first().copied().unwrap_or("?");
                let atoms_s = parts.get(1).copied().unwrap_or("?");
                let _ = writeln!(md, "| {} | {} | — | (×{}) |", sz, atoms_s, count);
            }
            let _ = writeln!(md);
        }
    }

    // ── Re-run ────────────────────────────────────────────────────────────────
    let _ = writeln!(md, "## Re-run (seed={})\n", rerun_seed);
    let form_b   = (initial_key >> 8) as u8;
    let _ = writeln!(md, "- **Initial atoms**: {} × (form={}, state=0)", n_atoms, form_b);
    let _ = writeln!(md, "- **tmax**: {}\n", tmax);

    let mut rng = StdRng::seed_from_u64(rerun_seed);
    let mut reactor = Reactor::new();
    for r in gene.to_reactions() { reactor.add_reaction(r); }
    reactor.add_atom_from_keys(initial_key, n_atoms);
    reactor.main_loop(tmax, &mut rng);

    let seq    = reactor.reaction_sequences().to_vec();
    let snaps  = reactor.state_snapshots().to_vec();
    let rstrs: Vec<String> = reactor.reactions().iter().map(|r| r.to_string()).collect();
    let mols   = reactor.find_molecules();
    let linked: Vec<&Molecule> = mols.iter().filter(|m| m.size() > 1).collect();

    let _ = writeln!(md, "**Reactions fired**: {}  |  **Simulation time reached**: {:.4}\n",
        seq.len(), reactor.current_time());

    // Reaction fire count summary
    if !seq.is_empty() {
        let mut fire_counts = vec![0u64; reactor.reactions().len()];
        for &(_, idx) in &seq { if idx < fire_counts.len() { fire_counts[idx] += 1; } }
        let _ = writeln!(md, "### Reaction Fire Counts\n");
        let _ = writeln!(md, "| Reaction | Count | % |");
        let _ = writeln!(md, "|----------|-------|---|");
        let total_f = seq.len() as f64;
        let mut fire_rows: Vec<(String, u64)> = rstrs.iter().cloned()
            .zip(fire_counts.iter().copied()).collect();
        fire_rows.sort_by(|a, b| b.1.cmp(&a.1));
        for (rstr, cnt) in &fire_rows {
            if *cnt > 0 {
                let _ = writeln!(md, "| `{}` | {} | {:.1}% |",
                    rstr, cnt, 100.0 * *cnt as f64 / total_f);
            }
        }
        let _ = writeln!(md);

        // Full sequence (first 200 steps) with per-step atom/link breakdown
        let show = seq.len().min(200);
        if seq.len() > 200 {
            let _ = writeln!(md,
                "*Showing first 200 of {} events — full sequence in `reaction_logs.jsonl`.*\n",
                seq.len());
        } else {
            let _ = writeln!(md, "### Reaction Sequence (re-run)\n");
        }
        let _ = writeln!(md, "| Step | Time | Reaction | Atoms (after) | Bonds (after) |");
        let _ = writeln!(md, "|------|------|----------|---------------|---------------|");
        for (i, &(t, idx)) in seq.iter().take(show).enumerate() {
            let (atoms_str, links_str) = snaps.get(i)
                .map(|(a, l)| format_snapshot(a, l))
                .unwrap_or_else(|| ("—".to_string(), "—".to_string()));
            let _ = writeln!(md, "| {} | {:.6} | `{}` | {} | {} |",
                i + 1, t,
                rstrs.get(idx).cloned().unwrap_or_default(),
                atoms_str, links_str);
        }
        let _ = writeln!(md);
    }

    // ── Final molecules ───────────────────────────────────────────────────────
    let _ = writeln!(md, "### Final Molecules (re-run)\n");
    let _ = writeln!(md, "| # | Size | Atoms (form,state) | Bonds | Comp. sim | Size ratio | Degree dist | WL dist (d=2) | Iso |");
    let _ = writeln!(md, "|---|------|--------------------|-------|-----------|------------|-------------|---------------|-----|");

    let mut mol_row = 0usize;
    for m in &mols {
        let mut sorted: Vec<(u8, u8)> = m.atoms.values().cloned().collect();
        sorted.sort();
        let atom_str = sorted.iter()
            .map(|&(f, s)| format!("({},{})", f, s))
            .collect::<Vec<_>>().join(" ");

        if m.size() > 1 {
            mol_row += 1;
            let comp  = m.composition_similarity(target);
            let sratio = m.size() as f64 / target.size().max(1) as f64;
            let deg   = m.degree_sequence_distance(target);
            let wl    = m.wl_distance(target, 2);
            let iso_s = if target.size() <= iso_max && m.size() == target.size() {
                if comp > 0.8 && m.is_isomorphic_to(target) { "✓" } else { "✗" }
            } else { "—" };
            let _ = writeln!(md, "| {} | {} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {} |",
                mol_row, m.size(), atom_str, m.links.len(),
                comp, sratio, deg, wl, iso_s);
        } else {
            let _ = writeln!(md, "| — | {} | {} | {} | — | — | — | — | — |",
                m.size(), atom_str, m.links.len());
        }
    }
    if linked.is_empty() {
        let _ = writeln!(md, "*No linked molecules produced.*");
    }
    let _ = writeln!(md);

    // ── Score breakdown ───────────────────────────────────────────────────────
    let bd = score_full(&reactor, target, iso_max);
    let _ = writeln!(md, "### Score Breakdown (re-run)\n");
    let _ = writeln!(md, "| Component | Raw value | Multiplier | Contribution |");
    let _ = writeln!(md, "|-----------|-----------|------------|-------------|");
    let _ = writeln!(md, "| Linked molecules | {} | × 0.5 | {:.4} |",
        bd.n_linked, bd.score_linked);
    let _ = writeln!(md, "| Best composition similarity | {:.4} | × 10 | {:.4} |",
        bd.best_comp, bd.score_comp);
    let _ = writeln!(md, "| Best size match | {:.4} | × 5 | {:.4} |",
        bd.best_size, bd.score_size);
    let iso_label = if bd.iso_found { "✓ exact match" } else { "✗ none" };
    let _ = writeln!(md, "| Isomorphic match | {} | +100 | {:.4} |",
        iso_label, bd.score_iso);
    let _ = writeln!(md, "| **TOTAL** | | | **{:.4}** |", bd.score_total);
    let _ = writeln!(md);

    // ── Target ────────────────────────────────────────────────────────────────
    let _ = writeln!(md, "## Target Molecule: `{}`\n", target.name);
    let _ = writeln!(md, "**Size**: {} atoms, {} bonds\n",
        target.atoms.len(), target.links.len());
    let _ = writeln!(md, "| Atom id | Form | State |");
    let _ = writeln!(md, "|---------|------|-------|");
    let mut tatoms: Vec<(u32, u8, u8)> = target.atoms.iter()
        .map(|(&id, &(f, s))| (id, f, s)).collect();
    tatoms.sort_by_key(|(id, _, _)| *id);
    for (id, f, s) in &tatoms {
        let _ = writeln!(md, "| {} | {} | {} |", id, f, s);
    }
    let _ = writeln!(md);
    let bond_strs: Vec<String> = target.links.values()
        .map(|(a, b)| format!("{a}—{b}")).collect();
    let _ = writeln!(md, "**Bonds**: {}\n",
        if bond_strs.is_empty() { "none".to_string() }
        else { bond_strs.join(", ") });

    // ── Similarity details for ALL linked molecules ───────────────────────────
    let _ = writeln!(md, "## Similarity Details — All Linked Molecules vs Target\n");
    if linked.is_empty() {
        let _ = writeln!(md, "*No linked molecules to compare.*\n");
    } else {
        let _ = writeln!(md, "| # | Size | Comp. similarity | Size ratio | Degree distance | WL distance (d=2) | Isomorphic |");
        let _ = writeln!(md, "|---|------|-----------------|------------|-----------------|-------------------|------------|");
        for (i, m) in linked.iter().enumerate() {
            let comp  = m.composition_similarity(target);
            let sratio = m.size() as f64 / target.size().max(1) as f64;
            let deg   = m.degree_sequence_distance(target);
            let wl    = m.wl_distance(target, 2);
            let iso_s = if target.size() <= iso_max && m.size() == target.size() {
                if comp > 0.8 && m.is_isomorphic_to(target) { "✓ yes" } else { "✗ no" }
            } else { "— (size mismatch)" };
            let _ = writeln!(md, "| {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {} |",
                i + 1, m.size(), comp, sratio, deg, wl, iso_s);
        }
        let _ = writeln!(md);
    }

    md
}

// ── Main ─────────────────────────────────────────────────────────────────────
fn main() {
    let cli = Cli::parse();

    // Load phylogeny
    let phylo_path = format!("{}/phylogeny.json", cli.data_dir);
    let phylo_str = std::fs::read_to_string(&phylo_path).unwrap_or_else(|e| {
        eprintln!("Cannot read {phylo_path}: {e}"); std::process::exit(1);
    });
    let phylo: serde_json::Value = serde_json::from_str(&phylo_str).unwrap_or_else(|e| {
        eprintln!("Cannot parse phylogeny.json: {e}"); std::process::exit(1);
    });

    // Find node
    let node = phylo["nodes"].as_array()
        .and_then(|arr| arr.iter().find(|n| n["id"].as_u64() == Some(cli.id as u64)))
        .unwrap_or_else(|| {
            eprintln!("Individual id={} not found in {phylo_path}", cli.id);
            std::process::exit(1);
        });

    // Config
    let cfg       = &phylo["config"];
    let max_count = cfg["max_count"].as_u64().unwrap_or(1) as u8;
    let n_atoms   = cfg["n_atoms"].as_u64().unwrap_or(10) as usize;
    let tmax      = cfg["tmax"].as_f64().unwrap_or(1000.0);
    let form_u8   = cfg["form"].as_u64().unwrap_or(1) as u8;
    let initial_key: u16 = (form_u8 as u16) << 8;
    let rerun_seed = if cli.seed == 0 {
        cfg["seed"].as_u64().unwrap_or(42)
    } else {
        cli.seed
    };

    // Target
    let target = build_target(); // target_from_json(&phylo["target"]);

    // Genome
    let genome_str = node["genome"].as_str().unwrap_or("");
    let genome: Vec<u8> = genome_str.chars()
        .map(|c| u8::from_str_radix(&c.to_string(), 16).unwrap_or(0))
        .collect();
    let gene = Gene8 { genome, max_count };

    // Original reaction log (optional)
    let log_path = format!("{}/reaction_logs.jsonl", cli.data_dir);
    let orig_log: Option<serde_json::Value> = std::fs::read_to_string(&log_path).ok()
        .and_then(|s| {
            s.lines().find_map(|line| {
                let v: serde_json::Value = serde_json::from_str(line).ok()?;
                if v["individual_id"].as_u64() == Some(cli.id as u64) { Some(v) } else { None }
            })
        });

    // Generate
    let md = generate_markdown(
        cli.id, node, &gene, &target,
        initial_key, n_atoms, tmax, 10,
        rerun_seed, orig_log.as_ref(),
    );

    // Write
    let out_path = cli.out.unwrap_or_else(||
        format!("{}/markdown/individual_{}.md", cli.data_dir, cli.id)
    );
    if let Some(p) = std::path::Path::new(&out_path).parent() {
        std::fs::create_dir_all(p).ok();
    }
    std::fs::write(&out_path, &md).unwrap_or_else(|e| {
        eprintln!("Cannot write {out_path}: {e}"); std::process::exit(1);
    });
    println!("Report written: {out_path}");
}
