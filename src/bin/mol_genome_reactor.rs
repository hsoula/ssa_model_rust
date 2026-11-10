//! mol_genome_reactor — spawn molecules from a MolGenome into a reactor and run it.
//!
//! The reactor is initialised with:
//!   • N copies of each molecule decoded from the genome (structurally inserted)
//!   • M free individual atoms of each form/state present in the decoded molecules
//!   • Reactions loaded from a chemistry JSON file (same format as gen_chemistry)
//!
//! Usage:
//!   mol_genome_reactor --genome E100201000F --chemistry data/chemistry.json
//!   mol_genome_reactor --random --size 120 --copies 3 --free-atoms 5 \
//!                      --chemistry data/chemistry.json --tmax 500 --seed 7

use clap::Parser;
use rand::{SeedableRng, rngs::StdRng};
use ssa_model::mol_genome::MolGenome;
use ssa_model::molecules::Molecule;
use ssa_model::reactions::Reaction;
use ssa_model::reactor::Reactor;
use std::collections::HashMap;

// ── CLI ───────────────────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(name = "mol_genome_reactor",
          about = "Spawn molecules from a MolGenome into a reactor and run a chemistry")]
struct Cli {
    /// Hex genome string (e.g. "E100201000F").  Omit to use --random.
    #[arg(long)]
    genome: Option<String>,

    /// Generate a random genome
    #[arg(long, default_value = "false")]
    random: bool,

    /// Length of random genome (tokens)
    #[arg(long, default_value = "120")]
    size: usize,

    /// Maximum atom form
    #[arg(long, default_value = "2")]
    max_form: u8,

    /// Maximum atom state
    #[arg(long, default_value = "4")]
    max_state: u8,

    /// Number of copies of each decoded molecule to insert
    #[arg(long, default_value = "3")]
    copies: usize,

    /// Number of free individual atoms per (form,state) pair found in the molecules
    #[arg(long, default_value = "5")]
    free_atoms: usize,

    /// Chemistry JSON file (from gen_chemistry)
    #[arg(long)]
    chemistry: Option<String>,

    /// Simulation end time
    #[arg(long, default_value = "500.0")]
    tmax: f64,

    /// Random seed
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Save final reactor state as markdown
    #[arg(long)]
    out_md: Option<String>,
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn atom_sym(form: u8, state: u8) -> String {
    format!("{}{}", (form + 96) as char, state)
}

fn rxn_type_label(r: &Reaction) -> &'static str {
    if r.is_production()  { return "production"; }
    if r.is_degradation() { return "degradation"; }
    match (r.is_substrate_contact(), r.is_product_contact()) {
        (true,  false) => "contact→bond",
        (true,  true)  => "contact",
        (false, true)  => "link→split",
        (false, false) => "link→state",
    }
}

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let cli = Cli::parse();
    let mut rng = StdRng::seed_from_u64(cli.seed);

    // ── 1. Build genome ───────────────────────────────────────────────────────
    let gene = if cli.random || cli.genome.is_none() {
        println!("Random genome  size={}  max_form={}  max_state={}  seed={}",
            cli.size, cli.max_form, cli.max_state, cli.seed);
        MolGenome::new(cli.size, cli.max_form, cli.max_state, &mut rng)
    } else {
        let s = cli.genome.unwrap();
        println!("Genome: {s}");
        MolGenome::from_hex(&s, cli.max_form, cli.max_state)
    };

    println!("  hex: {}", gene.to_hex());

    // ── 2. Decode molecules ───────────────────────────────────────────────────
    let mols = gene.to_molecules();
    println!("\n── Decoded molecules ({}) ─────────────────────────────────────", mols.len());
    if mols.is_empty() {
        println!("  (no valid molecules — at least 2 atoms required per pattern)");
        println!("  Hint: try --random or pass a longer --genome.");
    }
    for (i, m) in mols.iter().enumerate() {
        let mut atoms: Vec<(u32,u8,u8)> = m.atoms.iter()
            .map(|(&id, &(f,s))| (id,f,s)).collect();
        atoms.sort_by_key(|(id,_,_)| *id);
        let chain: Vec<String> = atoms.iter().map(|(_,f,s)| atom_sym(*f,*s)).collect();
        println!("  [{i}] {} atoms  {} bonds  chain: {}",
            m.atoms.len(), m.links.len(), chain.join("—"));
    }

    // ── 3. Collect (form,state) pairs present in decoded molecules ────────────
    let mut atom_types: HashMap<(u8,u8), usize> = HashMap::new();
    for m in &mols {
        for &(f, s) in m.atoms.values() {
            *atom_types.entry((f, s)).or_insert(0) += 1;
        }
    }

    // ── 4. Load chemistry ─────────────────────────────────────────────────────
    let reactions: Vec<Reaction> = match cli.chemistry.as_deref() {
        Some(path) => {
            let rxns = Reaction::news_from_json(path, 0);
            println!("\n── Chemistry: {} reactions from {path} ─────────────────────", rxns.len());
            rxns
        }
        None => {
            println!("\n── No chemistry file provided (--chemistry).  Running without reactions.");
            Vec::new()
        }
    };

    // ── 5. Build reactor ──────────────────────────────────────────────────────
    let mut reactor = Reactor::new();

    // Add reactions
    for r in &reactions {
        reactor.add_reaction(r.clone());
    }
    println!("\n── Reactor setup ─────────────────────────────────────────────────");

    // Insert molecule copies
    let mut mol_count = 0usize;
    for m in &mols {
        for _ in 0..cli.copies {
            reactor.add_molecule(m);
            mol_count += 1;
        }
    }
    println!("  Inserted {} molecule copies ({} distinct × {} copies)",
        mol_count, mols.len(), cli.copies);

    // Insert free individual atoms for each (form,state) present
    let mut free_count = 0usize;
    let mut sorted_types: Vec<(u8,u8)> = atom_types.keys().cloned().collect();
    sorted_types.sort();
    for (f, s) in &sorted_types {
        let key: u16 = (*f as u16) << 8 | (*s as u16);
        reactor.add_atom_from_keys(key, cli.free_atoms);
        free_count += cli.free_atoms;
        println!("  Free atoms: {} × {}{}",
            cli.free_atoms, (f + 96) as char, s);
    }
    println!("  Total free atoms: {free_count}");
    println!("  Total atoms at start: {}",
        mol_count * mols.first().map(|m| m.atoms.len()).unwrap_or(0) + free_count);

    // ── 6. Print reaction list ────────────────────────────────────────────────
    if !reactions.is_empty() {
        println!("\n── Reactions ─────────────────────────────────────────────────────");
        for r in &reactions {
            println!("  [{:>3}] {}  ({}  speed={:.3})",
                r.id(), r.to_string(), rxn_type_label(r), r.reaction_speed());
        }
    }

    // ── 7. Run simulation ─────────────────────────────────────────────────────
    println!("\n── Running simulation  tmax={} ────────────────────────────────────", cli.tmax);
    let t_end = reactor.main_loop(cli.tmax, &mut rng);
    let total_fired = reactor.reaction_sequences().len();
    println!("  Finished at t={t_end:.4}  reactions fired={total_fired}");

    // ── 8. Reaction fire counts ───────────────────────────────────────────────
    if total_fired > 0 && !reactions.is_empty() {
        println!("\n── Fire counts ──────────────────────────────────────────────────");
        println!("  {:>6}  {:>6}  {}  (type)", "fires", "%", "reaction");
        println!("  {:->6}  {:->6}  {:->40}", "", "", "");
        for (r, count) in reactor.reaction_fire_counts() {
            if count == 0 { continue; }
            let pct = 100.0 * count as f64 / total_fired as f64;
            println!("  {:>6}  {:>5.1}%  {}  ({})",
                count, pct, r.to_string(), rxn_type_label(r));
        }
    }

    // ── 9. Final atom distribution ────────────────────────────────────────────
    println!("\n── Final atom distribution ──────────────────────────────────────");
    let mut final_atoms: HashMap<u16, usize> = HashMap::new();
    for a in reactor.iter() {
        *final_atoms.entry(a.get_key()).or_insert(0) += 1;
    }
    let mut final_v: Vec<_> = final_atoms.iter().collect();
    final_v.sort_by_key(|(k, _)| **k);
    let total_atoms = reactor.get_atom_number();
    for (key, cnt) in &final_v {
        let f = ((*key >> 8) & 0xff) as u8;
        let s = (*key & 0xff) as u8;
        println!("  {}{}  ×{}  ({:.1}%)",
            (f + 96) as char, s, *cnt,
            100.0 * (**cnt) as f64 / total_atoms.max(1) as f64);
    }
    println!("  Total atoms: {total_atoms}");

    // ── 10. Final molecule distribution ───────────────────────────────────────
    let final_mols = reactor.find_molecules();
    let mut size_dist: HashMap<usize, usize> = HashMap::new();
    for m in &final_mols { *size_dist.entry(m.size()).or_insert(0) += 1; }
    let mut sizes: Vec<_> = size_dist.iter().collect();
    sizes.sort_by_key(|(s, _)| **s);

    println!("\n── Final molecule size distribution ({} components) ─────────────",
        final_mols.len());
    for (sz, cnt) in &sizes {
        let bar: String = "#".repeat((*cnt * 30 / final_mols.len().max(1)).max(1));
        println!("  size {:>3}  ×{:>4}  {}", sz, cnt, bar);
    }

    let linked: Vec<&Molecule> = final_mols.iter().filter(|m| m.size() > 1).collect();
    println!("  Linked molecules (size>1): {}", linked.len());

    if !linked.is_empty() {
        println!("\n── Linked molecules detail ───────────────────────────────────────");
        for (i, m) in linked.iter().enumerate() {
            let mut atoms: Vec<(u32,u8,u8)> = m.atoms.iter()
                .map(|(&id,&(f,s))| (id,f,s)).collect();
            atoms.sort_by_key(|(id,_,_)| *id);
            let chain: Vec<String> = atoms.iter()
                .map(|(_,f,s)| atom_sym(*f,*s)).collect();
            println!("  [{i:>3}] size={}  bonds={}  {}",
                m.atoms.len(), m.links.len(), chain.join("—"));
        }
    }

    // ── 11. Optional markdown dump ────────────────────────────────────────────
    if let Some(md_path) = cli.out_md {
        reactor.save_markdown(&md_path);
        println!("\nMarkdown saved → {md_path}");
    }
}
