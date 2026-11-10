//! simple_ga — single-form GA micro-experiment.
//!
//! One form (=1), states 0–2, target = dimer (1,1)—(1,2).
//! 10 generations × 10 individuals, 5 initial atoms of (form=1, state=0).
//!
//! Outputs (all in results/simple_ga/):
//!   ga_data.csv           — per-individual GA stats (all gens × all inds)
//!   reaction_logs.jsonl   — ordered reaction calls + final molecules per individual
//!   final_molecules.jsonl — molecules at end of each simulation
//!   phylogeny.json        — ancestor graph + point-mutation records

use rand::{SeedableRng, rngs::StdRng, Rng};
use serde::Serialize;
use ssa_model::Gene8::Gene8;
use ssa_model::molecules::Molecule;
use ssa_model::reactor::Reactor;
use std::collections::HashMap;
use std::io::Write;

// ── Experiment constants ──────────────────────────────────────────────────────
const FORM: u8            = 1;
const MAX_COUNT: u8       = 1;    // Gene8 form cap; atom states range 0..=MAX_COUNT*2 = 0..=2
const N_ATOMS: usize      = 10;   // 5 × target_size (target has 2 atoms) → 10 initial atoms
const GENOME_SIZE: usize  = 400;  // genome token length
const POPULATION: usize   = 20;
const GENERATIONS: usize  = 10;
const MUTATION_RATE: f64  = 0.1;
const TMAX: f64            = 1000.0;
const SEED: u64            = 42;
// One individual in the initial population is seeded with a known valid bonding reaction.
// With max_count=1, states 0-2, random genomes rarely contain a reaction that can fire from
// all-zero initial atoms; seeding provides a fitness gradient the GA can select on.
const N_SEED: usize       = 0;
const ISO_MAX_SIZE: usize  = 10;
const OUT: &str            = "results/simple_ga";

// ── Target: dimer (form=1,state=1) — (form=1,state=2) ────────────────────────
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

// ── Score breakdown ───────────────────────────────────────────────────────────
#[derive(Clone, Serialize)]
struct ScoreBreakdown {
    n_linked: usize,
    score_linked: f64,
    score_comp: f64,
    score_size: f64,
    score_iso: f64,
    score_total: f64,
}

fn score_detailed(reactor: &Reactor, target: &Molecule) -> ScoreBreakdown {
    let mols = reactor.find_molecules();
    let linked: Vec<&Molecule> = mols.iter().filter(|m| m.size() > 1).collect();
    let n_linked = linked.len();
    let t_size = target.size() as f64;

    let best_comp = linked.iter()
        .map(|m| m.composition_similarity(target))
        .fold(0.0_f64, f64::max);

    let best_size = if linked.is_empty() {
        0.0
    } else {
        linked.iter()
            .map(|m| 1.0 / (1.0 + (m.size() as f64 - t_size).abs() / t_size.max(1.0)))
            .fold(0.0_f64, f64::max)
    };

    let iso = if target.size() <= ISO_MAX_SIZE {
        linked.iter()
            .filter(|m| m.size() == target.size())
            .any(|m| m.composition_similarity(target) > 0.8 && m.is_isomorphic_to(target))
    } else {
        false
    };

    let score_linked = n_linked as f64 / 2.0;
    let score_comp   = best_comp * 10.0;
    let score_size   = best_size * 5.0;
    let score_iso    = if iso { 100.0 } else { 0.0 };

    ScoreBreakdown {
        n_linked,
        score_linked,
        score_comp,
        score_size,
        score_iso,
        score_total: score_linked + score_comp + score_size + score_iso,
    }
}

// ── Serialisable records ──────────────────────────────────────────────────────
#[derive(Serialize)]
struct ReactionEvent {
    time: f64,
    reaction_idx: usize,
    reaction_str: String,
    /// Atom counts after the reaction: symbol → count  (e.g. {"a0": 8, "a1": 1})
    atoms: HashMap<String, usize>,
    /// Bond counts after the reaction (normalised key): symbol → count (e.g. {"a1-a2": 1})
    links: HashMap<String, usize>,
}

#[derive(Clone, Serialize)]
struct MoleculeRecord {
    size: usize,
    atoms: Vec<(u8, u8)>,  // (form, state) sorted
    n_links: usize,
}

fn mol_record(m: &Molecule) -> MoleculeRecord {
    let mut atoms: Vec<(u8, u8)> = m.atoms.values().cloned().collect();
    atoms.sort();
    MoleculeRecord { size: m.size(), atoms, n_links: m.links.len() }
}

#[derive(Serialize)]
struct SimLog {
    individual_id: usize,
    generation: usize,
    position: usize,
    total_reactions: usize,
    reaction_sequence: Vec<ReactionEvent>,
    final_molecules: Vec<MoleculeRecord>,
}

// ── GA helpers ────────────────────────────────────────────────────────────────
fn tournament_idx(scores: &[f64], k: usize, rng: &mut StdRng) -> usize {
    let mut best = rng.random_range(0..scores.len());
    for _ in 1..k {
        let idx = rng.random_range(0..scores.len());
        if scores[idx] > scores[best] { best = idx; }
    }
    best
}

/// Positions where two equal-length byte slices differ: (pos, old, new).
fn point_diff(before: &[u8], after: &[u8]) -> Vec<(usize, u8, u8)> {
    before.iter().zip(after.iter()).enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, (&a, &b))| (i, a, b))
        .collect()
}

// ── Snapshot helpers ──────────────────────────────────────────────────────────
/// Convert atom key (form<<8 | state) to symbol like "a0", "a1", "b2".
fn key_to_sym(k: u16) -> String {
    let form  = ((k >> 8) & 0xff) as u8;
    let state = (k & 0xff) as u8;
    format!("{}{}", (form + 96) as char, state)
}

/// Normalised link key as "aX-aY" (lexicographically smaller first).
fn normalize_link_key(k: u32) -> String {
    let a = key_to_sym((k >> 16) as u16);
    let b = key_to_sym((k & 0xffff) as u16);
    if a <= b { format!("{a}-{b}") } else { format!("{b}-{a}") }
}

// ── Seed individual ───────────────────────────────────────────────────────────
/// Return a Gene8 with a known-valid bonding reaction embedded at the start.
///
/// Embeds: a0+a0 → a1.a2  (substrate-contact, creates bond, changes states 0→1 and 0→2)
/// Gene8 encoding:
///   E(0E) sf1++(00) sf2++(04) contact++(08) contact++(08) ps1++(09) ps2++(0B) ps2++(0C) F(0F)
/// With max_count=1 this decodes to:
///   substrate: (form=1,state=0) + (form=1,state=0)  →  product: (form=1,state=1) . (form=1,state=2)
fn seed_genome(size: usize, rng: &mut StdRng) -> Gene8 {
    let seed_pattern: &[u8] = &[0x0E, 0x00, 0x04, 0x08, 0x08, 0x09, 0x0B, 0x0C, 0x0F];
    let mut g = Gene8::new_with_max(size, MAX_COUNT, rng);
    for (i, &v) in seed_pattern.iter().enumerate() {
        g.genome[i] = v;
    }
    g
}

// ── Main ─────────────────────────────────────────────────────────────────────
fn main() {
    let mut rng = StdRng::seed_from_u64(SEED);
    std::fs::create_dir_all(OUT).unwrap();

    let target = build_target();
    let initial_key: u16 = (FORM as u16) << 8; // (form=1, state=0)

    println!("Target: {} atoms, {} bonds  [form={}, states 0..{}]",
        target.size(), target.links.len(), FORM, MAX_COUNT * 2);
    println!("Initial atoms: {} × key=0x{:04X}  tmax={}  seed={}",
        N_ATOMS, initial_key, TMAX, SEED);
    println!("Population={}  Generations={}  genome_size={}\n",
        POPULATION, GENERATIONS, GENOME_SIZE);

    // ── Open output files ─────────────────────────────────────────────────────
    let mut csv_file = std::fs::File::create(format!("{OUT}/ga_data.csv")).unwrap();
    writeln!(csv_file,
        "generation,individual_id,position,parent1_id,parent2_id,\
         genome_len,n_reactions_decoded,total_reactions_called,\
         score_total,score_linked,score_comp,score_size,score_iso,\
         n_linked_molecules,is_elite,n_point_mutations"
    ).unwrap();

    let mut rxn_file = std::fs::File::create(format!("{OUT}/reaction_logs.jsonl")).unwrap();
    let mut mol_file = std::fs::File::create(format!("{OUT}/final_molecules.jsonl")).unwrap();

    // Phylogeny accumulators
    let mut phyl_nodes: Vec<serde_json::Value> = Vec::new();
    let mut phyl_edges: Vec<serde_json::Value> = Vec::new();
    let mut phyl_muts:  Vec<serde_json::Value> = Vec::new();

    // ── Initial population ────────────────────────────────────────────────────
    // The first N_SEED individuals carry a known valid bonding reaction;
    // the remainder are fully random.  This provides a fitness gradient from
    // generation 0 so the GA has something to select on.
    let mut pop: Vec<Gene8> = (0..POPULATION)
        .map(|i| if i < N_SEED {
            seed_genome(GENOME_SIZE, &mut rng)
        } else {
            Gene8::new_with_max(GENOME_SIZE, MAX_COUNT, &mut rng)
        })
        .collect();

    // Per-slot metadata: (individual_id, parent1_id, parent2_id, genome_after_crossover, is_elite)
    // For initial pop: genome_after_crossover == final genome, no parents.
    let mut pop_meta: Vec<(usize, Option<usize>, Option<usize>, String, bool)> =
        (0..POPULATION).map(|i| (i, None, None, pop[i].to_string(), false)).collect();

    let mut next_id: usize = POPULATION;
    let mut best_ever_score = f64::NEG_INFINITY;

    // ── GA loop ───────────────────────────────────────────────────────────────
    for generation in 0..GENERATIONS {
        let mut scores: Vec<f64> = Vec::new();

        // ── Evaluate all individuals ──────────────────────────────────────────
        for pos in 0..POPULATION {
            let (ind_id, p1_id, p2_id, ref genome_raw, is_elite) = pop_meta[pos].clone();
            let genome_raw = genome_raw.clone();

            // Build and run reactor
            let mut reactor = Reactor::new();
            for r in pop[pos].to_reactions() { reactor.add_reaction(r); }
            reactor.add_atom_from_keys(initial_key, N_ATOMS);
            reactor.main_loop(TMAX, &mut rng);

            let bd = score_detailed(&reactor, &target);
            let n_decoded = pop[pos].to_reactions().len();
            let total_calls = reactor.reaction_sequences().len();
            scores.push(bd.score_total);

            // n_point_mutations = positions that differ between genome_raw and final genome
            let genome_final = pop[pos].to_string();
            let n_pmut = {
                let rb: Vec<u8> = genome_raw.chars()
                    .map(|c| u8::from_str_radix(&c.to_string(), 16).unwrap_or(0))
                    .collect();
                let rf: Vec<u8> = genome_final.chars()
                    .map(|c| u8::from_str_radix(&c.to_string(), 16).unwrap_or(0))
                    .collect();
                if rb.len() == rf.len() { point_diff(&rb, &rf).len() } else { 0 }
            };

            // CSV row
            writeln!(csv_file,
                "{},{},{},{},{},{},{},{},{:.4},{:.4},{:.4},{:.4},{:.4},{},{},{}",
                generation, ind_id, pos,
                p1_id.map(|x| x.to_string()).unwrap_or_default(),
                p2_id.map(|x| x.to_string()).unwrap_or_default(),
                pop[pos].len(), n_decoded, total_calls,
                bd.score_total, bd.score_linked, bd.score_comp, bd.score_size, bd.score_iso,
                bd.n_linked, is_elite as u8, n_pmut
            ).unwrap();

            // Reaction sequence log (with per-step atom/link snapshot)
            let rxn_strs: Vec<String> = reactor.reactions().iter()
                .map(|r| r.to_string()).collect();
            let snaps = reactor.state_snapshots().to_vec();
            let seq: Vec<ReactionEvent> = reactor.reaction_sequences().iter()
                .enumerate()
                .map(|(i, &(t, idx))| {
                    let (atom_map, link_map) = snaps.get(i).cloned().unwrap_or_default();
                    let atoms: HashMap<String, usize> = atom_map.iter()
                        .map(|(&k, &v)| (key_to_sym(k), v))
                        .collect();
                    let mut links: HashMap<String, usize> = HashMap::new();
                    for (&k, &v) in &link_map {
                        *links.entry(normalize_link_key(k)).or_insert(0) += v;
                    }
                    ReactionEvent {
                        time: t,
                        reaction_idx: idx,
                        reaction_str: rxn_strs.get(idx).cloned().unwrap_or_default(),
                        atoms,
                        links,
                    }
                }).collect();
            let final_mols: Vec<MoleculeRecord> = reactor.find_molecules().iter()
                .map(|m| mol_record(m)).collect();

            writeln!(rxn_file, "{}", serde_json::to_string(&SimLog {
                individual_id: ind_id,
                generation, position: pos,
                total_reactions: total_calls,
                reaction_sequence: seq,
                final_molecules: final_mols.clone(),
            }).unwrap()).unwrap();

            // Molecule log (compact)
            writeln!(mol_file, "{}", serde_json::json!({
                "individual_id": ind_id,
                "generation": generation,
                "position": pos,
                "score_total": bd.score_total,
                "molecules": serde_json::to_value(&final_mols).unwrap(),
            })).unwrap();

            // Phylogeny node
            phyl_nodes.push(serde_json::json!({
                "id": ind_id,
                "generation": generation,
                "position": pos,
                "genome": genome_final,
                "genome_after_crossover": genome_raw,
                "score": bd.score_total,
                "n_decoded_reactions": n_decoded,
                "is_elite": is_elite,
                "parent1_id": p1_id,
                "parent2_id": p2_id,
            }));

            if bd.score_total > best_ever_score { best_ever_score = bd.score_total; }

            println!(
                "gen {:2}  pos {:2}  id={:3}  score={:7.2}  rxn_dec={:2}  calls={:5}  linked={:2}",
                generation, pos, ind_id, bd.score_total, n_decoded, total_calls, bd.n_linked
            );
        }

        let best_idx = scores.iter().enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i).unwrap();

        println!("--- gen {:2}  best={:.2}  mean={:.2}\n",
            generation,
            scores[best_idx],
            scores.iter().sum::<f64>() / scores.len() as f64
        );

        // ── Breed next generation ─────────────────────────────────────────────
        let mut next_pop:  Vec<Gene8> = Vec::new();
        let mut next_meta: Vec<(usize, Option<usize>, Option<usize>, String, bool)> = Vec::new();

        // Elite: carry best individual forward with a new ID
        let elite_parent_id = pop_meta[best_idx].0;
        let elite_id = next_id; next_id += 1;
        next_pop.push(pop[best_idx].clone());
        next_meta.push((elite_id, Some(elite_parent_id), None, pop[best_idx].to_string(), true));
        phyl_edges.push(serde_json::json!({
            "child": elite_id, "parent": elite_parent_id, "role": "elite"
        }));

        // Fill remainder via crossover + point mutation
        while next_pop.len() < POPULATION {
            let p1i = tournament_idx(&scores, 3, &mut rng);
            let p2i = tournament_idx(&scores, 3, &mut rng);
            let p1_id = pop_meta[p1i].0;
            let p2_id = pop_meta[p2i].0;

            let (raw1, raw2) = pop[p1i].crossing(&pop[p2i], &mut rng);
            let child1 = raw1.point_mutations(MUTATION_RATE, &mut rng);

            let id1 = next_id; next_id += 1;
            for (mpos, old, nv) in point_diff(&raw1.genome, &child1.genome) {
                phyl_muts.push(serde_json::json!({
                    "individual_id": id1, "position": mpos, "old": old, "new": nv
                }));
            }
            phyl_edges.push(serde_json::json!({
                "child": id1, "parent1": p1_id, "parent2": p2_id, "role": "crossover"
            }));
            next_pop.push(child1);
            next_meta.push((id1, Some(p1_id), Some(p2_id), raw1.to_string(), false));

            if next_pop.len() < POPULATION {
                let child2 = raw2.point_mutations(MUTATION_RATE, &mut rng);
                let id2 = next_id; next_id += 1;
                for (mpos, old, nv) in point_diff(&raw2.genome, &child2.genome) {
                    phyl_muts.push(serde_json::json!({
                        "individual_id": id2, "position": mpos, "old": old, "new": nv
                    }));
                }
                phyl_edges.push(serde_json::json!({
                    "child": id2, "parent1": p1_id, "parent2": p2_id, "role": "crossover"
                }));
                next_pop.push(child2);
                next_meta.push((id2, Some(p1_id), Some(p2_id), raw2.to_string(), false));
            }
        }

        pop      = next_pop;
        pop_meta = next_meta;
    }

    // ── Save phylogeny ────────────────────────────────────────────────────────
    let phylogeny = serde_json::json!({
        "config": {
            "form": FORM, "max_count": MAX_COUNT, "n_atoms": N_ATOMS,
            "genome_size": GENOME_SIZE, "population": POPULATION,
            "generations": GENERATIONS, "mutation_rate": MUTATION_RATE,
            "tmax": TMAX, "seed": SEED,
        },
        "target": {
            "name": target.name,
            "atoms": target.atoms.iter().map(|(&id, &(f, s))| serde_json::json!({
                "id": id, "form": f, "state": s
            })).collect::<Vec<_>>(),
            "links": target.links.values().map(|(a, b)| serde_json::json!({
                "a": a, "b": b
            })).collect::<Vec<_>>(),
        },
        "nodes": phyl_nodes,
        "edges": phyl_edges,    
        "mutations": phyl_muts,
    });
    std::fs::write(
        format!("{OUT}/phylogeny.json"),
        serde_json::to_string_pretty(&phylogeny).unwrap(),
    ).unwrap();

    println!("=== Done  best_score_ever={:.2} ===", best_ever_score);
    println!("Output: {OUT}/");
    println!("  ga_data.csv           per-individual GA stats");
    println!("  reaction_logs.jsonl   ordered reaction calls + final molecules");
    println!("  final_molecules.jsonl molecules at simulation end");
    println!("  phylogeny.json        ancestor graph + point-mutation records");
}
