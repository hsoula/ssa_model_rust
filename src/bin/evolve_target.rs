use clap::Parser;
use rand::{SeedableRng, rngs::StdRng, Rng};
use ssa_model::Gene8::Gene8;
use ssa_model::molecules::Molecule;
use ssa_model::reactions::Reaction;
use ssa_model::reactor::Reactor;
use std::collections::HashSet;
use std::io::Write;
use std::path::PathBuf;

// ── CLI ──────────────────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "evolve_target",
    about = "Evolve a chemistry to produce target molecules from unlinked atoms in state 0"
)]
struct Cli {
    /// Target molecule JSON files (repeat: -m mol1.json -m mol2.json)
    #[arg(short = 'm', long = "molecule", required = true, num_args = 1..)]
    molecules: Vec<PathBuf>,

    /// Number of atoms per type to start with (all unlinked, state 0)
    #[arg(long, default_value = "100")]
    n_atoms: usize,

    /// Initial genome size (Gene8 tokens)
    #[arg(long, default_value = "200")]
    genome_size: usize,

    /// Population size
    #[arg(long, default_value = "20")]
    population: usize,

    /// Number of GA generations
    #[arg(long, default_value = "100")]
    generations: usize,

    /// Point mutation probability per genome token
    #[arg(long, default_value = "0.02")]
    mutation_rate: f64,

    /// Simulation time per individual
    #[arg(long, default_value = "2000")]
    tmax: f64,

    /// Max target molecule size to attempt graph isomorphism check (VF2 safety cap)
    #[arg(long, default_value = "20")]
    iso_max_size: usize,

    /// Random seed
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Output directory
    #[arg(short, long, default_value = "data/evolve_target")]
    output_dir: PathBuf,
}

// ── Atom initialisation ───────────────────────────────────────────────────────

/// Remap reaction atom forms onto the actual target atom forms.
///
/// Gene8 produces form values 1..max_count where max_count = forms.len().
/// This maps form value F → forms[F-1], a direct 1:1 correspondence.
/// Form 0 is left as-is (invalid, will be filtered by to_reactions).
fn remap_reactions(reactions: Vec<Reaction>, forms: &[u8]) -> Vec<Reaction> {
    if forms.is_empty() { return reactions; }

    let remap_key = |key: u16| -> u16 {
        let form  = (key >> 8) as u8;
        let state = (key & 0xFF) as u8;
        if form == 0 { return key; }
        let new_form = forms[((form as usize).wrapping_sub(1)) % forms.len()];
        ((new_form as u16) << 8) | (state as u16)
    };

    reactions.into_iter().map(|r| {
        let s1 = remap_key(r.get_substrate_a());
        let s2 = remap_key(r.get_substrate_b());
        let p1 = remap_key(r.get_product_a());
        let p2 = remap_key(r.get_product_b());
        let sk = ((s1 as u32) << 16) | (s2 as u32);
        let pk = ((p1 as u32) << 16) | (p2 as u32);
        r.set_substrate_key(sk).set_product_key(pk)
    }).collect()
}

/// Collect every distinct atom form (type) appearing in the target molecules.
fn forms_from_targets(targets: &[Molecule]) -> Vec<u8> {
    let mut seen: HashSet<u8> = HashSet::new();
    for mol in targets {
        for &(form, _) in mol.atoms.values() {
            if form > 0 { seen.insert(form); }
        }
    }
    let mut v: Vec<u8> = seen.into_iter().collect();
    v.sort();
    v
}

// ── Fitness ───────────────────────────────────────────────────────────────────

/// Score a reactor's connected components against the target list.
///
/// Per target:
///   +10  × best composition_similarity found in any linked molecule
///   + 5  × 1/(1 + relative_size_error) of the closest-size molecule
///   +100   if an exact isomorphic match exists (capped by iso_max_size)
fn score_reactor(reactor: &Reactor, targets: &[Molecule], iso_max_size: usize) -> f64 {
    let mols = reactor.find_molecules();
    let linked: Vec<&Molecule> = mols.iter().filter(|m| m.size() > 1).collect();

    let mut total = (linked.len() as f64) / 2.0;
    for target in targets {
        if linked.is_empty() { continue; }

        let t_size = target.size() as f64;

        // Composition similarity (0..1)
        let best_comp: f64 = linked.iter()
            .map(|m| m.composition_similarity(target))
            .fold(0.0_f64, f64::max);
        total += best_comp * 10.0;

        // Size proximity (0..1)
        let best_size: f64 = linked.iter()
            .map(|m| 1.0 / (1.0 + (m.size() as f64 - t_size).abs() / t_size.max(1.0)))
            .fold(0.0_f64, f64::max);
        total += best_size * 5.0;

        // Exact topology (expensive – guarded by size cap and composition threshold)
        if target.size() <= iso_max_size {
            for m in linked.iter().filter(|m| m.size() == target.size()) {
                if m.composition_similarity(target) > 0.8 && m.is_isomorphic_to(target) {
                    total += 100.0;
                    break;
                }
            }
        }
    }
    total
}

/// Build a fresh reactor from a Gene8 genome and evaluate it.
fn evaluate(
    gene:         &Gene8,
    initial_keys: &[u16],
    forms:        &[u8],
    n_atoms:      usize,
    targets:      &[Molecule],
    tmax:         f64,
    iso_max_size: usize,
    rng:          &mut StdRng,
) -> f64 {
    let raw = gene.to_reactions();
    if raw.is_empty() { return 0.0; }
    let reactions = remap_reactions(raw, forms);

    let mut reactor = Reactor::new();
    for r in reactions { reactor.add_reaction(r); }
    for &key in initial_keys { reactor.add_atom_from_keys(key, n_atoms); }
    reactor.main_loop(tmax, rng);

    score_reactor(&reactor, targets, iso_max_size)
}

// ── GA helpers ────────────────────────────────────────────────────────────────

fn tournament<'a>(pop: &'a [Gene8], scores: &[f64], k: usize, rng: &mut StdRng) -> &'a Gene8 {
    let mut best = rng.random_range(0..pop.len());
    for _ in 1..k {
        let idx = rng.random_range(0..pop.len());
        if scores[idx] > scores[best] { best = idx; }
    }
    &pop[best]
}
    
// ── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let cli = Cli::parse();
    let mut rng = StdRng::seed_from_u64(cli.seed);
    std::fs::create_dir_all(&cli.output_dir).unwrap();

    // Load targets
    let targets: Vec<Molecule> = cli.molecules.iter()
        .map(|p| {
            Molecule::load_from_json_file(p.to_str().unwrap()).unwrap_or_else(|e| {
                eprintln!("Cannot load {}: {e}", p.display());
                std::process::exit(1);
            })
        })
        .collect();

    println!("=== Target molecules ===");
    for mol in &targets {
        println!("  {} ({})  {} atoms  {} bonds",
            mol.name, mol.formula, mol.size(), mol.links.len());
    }

    // Initial atom pool: one key per form, all in state 0
    let forms = forms_from_targets(&targets);
    let initial_keys: Vec<u16> = forms.iter().map(|&f| (f as u16) << 8).collect();

    println!("\nAtom pool (state=0): {}",
        initial_keys.iter().map(|k| format!("0x{:04X}", k)).collect::<Vec<_>>().join(", "));
    println!("Count per type: {}  |  tmax: {}  |  iso_max_size: {}",
        cli.n_atoms, cli.tmax, cli.iso_max_size);
    println!("Population: {}  |  Generations: {}  |  Genome size: {}\n",
        cli.population, cli.generations, cli.genome_size);

    // Open fitness log
    let log_path = cli.output_dir.join("fitness_log.csv");
    let mut log = std::fs::File::create(&log_path).unwrap();
    writeln!(log, "generation,best,mean,worst,n_reactions,genome_len").unwrap();

    // Open JSON history (written at end)
    let mut history: Vec<serde_json::Value> = Vec::new();

    // max_count = number of distinct target forms so Gene8 form values map 1:1 onto them
    let max_count = forms.len() as u8;
    println!("Gene8 max_count: {} (distinct target forms: {:?})", max_count, forms);

    // Initial population
    let mut pop: Vec<Gene8> = (0..cli.population)
        .map(|_| Gene8::new_with_max(cli.genome_size, max_count, &mut rng))
        .collect();

    let mut best_ever_score = f64::NEG_INFINITY;
    let mut best_ever: Option<Gene8> = None;

    for generation in 0..cli.generations {
        // Evaluate
        let scores: Vec<f64> = pop.iter()
            .map(|g| evaluate(g, &initial_keys, &forms, cli.n_atoms, &targets,
                              cli.tmax, cli.iso_max_size, &mut rng))
            .collect();

        let best_idx = scores.iter().enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i).unwrap();

        let best_score  = scores[best_idx];
        let mean_score  = scores.iter().sum::<f64>() / scores.len() as f64;
        let worst_score = scores.iter().cloned().fold(f64::INFINITY, f64::min);
        let n_reactions = pop[best_idx].to_reactions().len();
        let genome_len  = pop[best_idx].len();

        println!(
            "gen {:3}  best={:7.2}  mean={:7.2}  worst={:7.2}  rxn={:3}  len={}",
            generation, best_score, mean_score, worst_score, n_reactions, genome_len
        );

        // CSV row
        writeln!(log, "{},{:.4},{:.4},{:.4},{},{}",
            generation, best_score, mean_score, worst_score, n_reactions, genome_len).unwrap();

        // JSON history entry
        history.push(serde_json::json!({
            "generation": generation,
            "best":       best_score,
            "mean":       mean_score,
            "worst":      worst_score,
            "n_reactions": n_reactions,
            "genome_len": genome_len,
        }));

        // Save new best
        if best_score > best_ever_score {
            best_ever_score = best_score;
            best_ever = Some(Gene8 { genome: pop[best_idx].genome.clone(), max_count: pop[best_idx].max_count });

            std::fs::write(
                cli.output_dir.join(format!("best_gen{:03}.genome", generation)),
                pop[best_idx].to_string(),
            ).unwrap();

            pop[best_idx].save_to_json(
                cli.output_dir.join("best_reactions.json").to_str().unwrap()
            );
            println!("  ↑ new best! ({:.2})", best_score);
        }

        // Breed next generation (elitism: carry best forward)
        let elite = Gene8 { genome: pop[best_idx].genome.clone(), max_count: pop[best_idx].max_count };
        let mut next = vec![elite];
        while next.len() < cli.population {
            let p1 = tournament(&pop, &scores, 3, &mut rng);
            let p2 = tournament(&pop, &scores, 3, &mut rng);
            let (c1, c2) = p1.crossing(p2, &mut rng);
            next.push(c1.point_mutations(cli.mutation_rate, &mut rng));
            if next.len() < cli.population {
                next.push(c2.point_mutations(cli.mutation_rate, &mut rng));
            }
        }
        pop = next;
    }

    // Save JSON fitness history
    let hist_path = cli.output_dir.join("fitness_history.json");
    std::fs::write(&hist_path, serde_json::to_string_pretty(&history).unwrap()).unwrap();

    // Final evaluation of the best individual
    if let Some(best) = &best_ever {
        let final_json = cli.output_dir.join("best_final.json");
        best.save_to_json(final_json.to_str().unwrap());

        // Longer final run (2× tmax) for the report
        let mut rng2 = StdRng::seed_from_u64(cli.seed);
        let mut reactor = Reactor::new();
        for r in best.to_reactions() { reactor.add_reaction(r); }
        for &key in &initial_keys { reactor.add_atom_from_keys(key, cli.n_atoms); }
        reactor.main_loop(cli.tmax * 2.0, &mut rng2);

        let final_score = score_reactor(&reactor, &targets, cli.iso_max_size);

        // Markdown report
        let md_path = cli.output_dir.join("best_final_report.md");
        let mut md = reactor.to_markdown();

        // Append target comparison section
        md.push_str("\n## Target comparison\n\n");
        let mols = reactor.find_molecules();
        for target in &targets {
            md.push_str(&format!("### {} ({})\n\n", target.name, target.formula));
            md.push_str(&format!("Target: {} atoms, {} bonds\n\n", target.size(), target.links.len()));

            let linked: Vec<&Molecule> = mols.iter().filter(|m| m.size() > 1).collect();
            if linked.is_empty() {
                md.push_str("*No linked molecules produced.*\n\n");
                continue;
            }

            let mut best_mols: Vec<(&Molecule, f64)> = linked.iter()
                .map(|m| (*m, m.composition_similarity(target)))
                .collect();
            best_mols.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

            md.push_str("| Rank | Size | Comp. similarity | Isomorphic? |\n");
            md.push_str("|------|------|-----------------|-------------|\n");
            for (rank, (mol, sim)) in best_mols.iter().take(5).enumerate() {
                let iso = if target.size() <= cli.iso_max_size && mol.size() == target.size() {
                    if mol.is_isomorphic_to(target) { "✓ yes" } else { "no" }
                } else { "—" };
                md.push_str(&format!(
                    "| {:4} | {:4} | {:>15.3} | {} |\n",
                    rank + 1, mol.size(), sim, iso
                ));
            }
            md.push('\n');
        }

        std::fs::write(&md_path, &md).unwrap();

        println!("\n=== Evolution complete ===");
        println!("Best score : {:.2}  (final re-run: {:.2})", best_ever_score, final_score);
        println!("Chemistry  : {}", final_json.display());
        println!("Report     : {}", md_path.display());
        println!("Fitness log: {}", log_path.display());
        println!("History    : {}", hist_path.display());
        println!("\nTo re-run the best chemistry:");
        let atom_args: Vec<String> = initial_keys.iter()
            .map(|k| format!("-a 0x{:04X}:{}", k, cli.n_atoms))
            .collect();
        println!(
            "  cargo run --release --bin simulate -- -r {} --tmax {} {}",
            final_json.display(),
            cli.tmax * 2.0,
            atom_args.join(" ")
        );
    }
}
