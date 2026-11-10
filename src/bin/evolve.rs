use clap::Parser;
use rand::{SeedableRng, rngs::StdRng, Rng};
use ssa_model::Gene8::Gene8;
use ssa_model::reactor::Reactor;
use ssa_model::utils::parse_atom_spec;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "evolve", about = "Evolve a reaction chemistry via genetic algorithm")]
struct Cli {
    /// Initial atoms as KEY:COUNT (e.g. 0x0100:300); repeat for multiple types.
    /// Either byte of the key may contain a '*' to wildcard it, e.g. 0x0*0*:50 (every
    /// form × every state), 0x0*03:50 (every form, state=3), 0x01**:50 (form=1, every state).
    #[arg(short, long = "atoms", required = true, num_args = 1..)]
    atoms: Vec<String>,

    /// Max form used to expand a wildcarded (*) form in an atom spec (forms 1..=max_form)
    #[arg(long, default_value = "10")]
    max_form: u8,

    /// Max state used to expand a wildcarded (*) state in an atom spec (states 0..=max_state)
    #[arg(long, default_value = "10")]
    max_state: u8,

    /// Genome size (number of Gene8 tokens per individual)
    #[arg(long, default_value = "200")]
    genome_size: usize,

    /// Population size
    #[arg(long, default_value = "20")]
    population: usize,

    /// Number of generations
    #[arg(long, default_value = "50")]
    generations: usize,

    /// Point mutation probability per token
    #[arg(long, default_value = "0.02")]
    mutation_rate: f64,

    /// Max simulation time per individual (keep short for GA speed)
    #[arg(long, default_value = "1000")]
    tmax: f64,

    /// Random seed
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Output directory (best genome + reactions JSON saved each generation)
    #[arg(short, long, default_value = "data/evolve")]
    output_dir: PathBuf,
}

/// Build and run a reactor from a Gene8 genome, return fitness score.
fn evaluate(gene: &Gene8, atom_specs: &[(u16, usize)], tmax: f64, rng: &mut StdRng) -> f64 {
    let reactions = gene.to_reactions();
    if reactions.is_empty() {
        return 0.0;
    }

    let mut reactor = Reactor::new();
    for r in reactions {
        reactor.add_reaction(r);
    }
    for &(key, count) in atom_specs {
        reactor.add_atom_from_keys(key, count);
    }

    reactor.main_loop(tmax, rng);

    let molecules = reactor.find_molecules();
    if molecules.is_empty() {
        return 0.0;
    }

    // Fitness: reward linked atoms and diversity of molecule sizes
    let n_linked: usize = molecules.iter()
        .filter(|m| m.size() > 1)
        .map(|m| m.size())
        .sum();
    let size_dist: HashMap<usize, usize> = molecules.iter()
        .map(|m| m.size())
        .fold(HashMap::new(), |mut acc, s| { *acc.entry(s).or_insert(0) += 1; acc });
    let n_distinct = size_dist.len() as f64;
    let max_size = *size_dist.keys().max().unwrap_or(&0) as f64;

    n_linked as f64 + n_distinct * 5.0 + max_size * 2.0
}

/// Tournament selection: pick best of k random individuals.
fn tournament<'a>(pop: &'a [Gene8], scores: &[f64], k: usize, rng: &mut StdRng) -> &'a Gene8 {
    let n = pop.len();
    let mut best_idx = rng.random_range(0..n);
    for _ in 1..k {
        let idx = rng.random_range(0..n);
        if scores[idx] > scores[best_idx] {
            best_idx = idx;
        }
    }
    &pop[best_idx]
}

fn main() {
    let cli = Cli::parse();
    let mut rng = StdRng::seed_from_u64(cli.seed);

    let atom_specs: Vec<(u16, usize)> = cli.atoms.iter()
        .flat_map(|s| parse_atom_spec(s, cli.max_form, cli.max_state).unwrap_or_else(|e| {
            eprintln!("Error parsing atom spec '{s}': {e}");
            std::process::exit(1);
        }))
        .collect();
    std::fs::create_dir_all(&cli.output_dir).unwrap();

    // Initial population
    let mut population: Vec<Gene8> = (0..cli.population)
        .map(|_| Gene8::new(cli.genome_size, &mut rng))
        .collect();

    let mut best_ever_score = f64::NEG_INFINITY;
    let mut best_ever_genome: Option<Gene8> = None;

    for generation in 0..cli.generations {
        // Evaluate
        let scores: Vec<f64> = population.iter()
            .map(|g| evaluate(g, &atom_specs, cli.tmax, &mut rng))
            .collect();

        let best_idx = scores.iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        let worst_score = scores.iter().cloned().fold(f64::INFINITY, f64::min);
        let mean_score = scores.iter().sum::<f64>() / scores.len() as f64;
        let best_score = scores[best_idx];
        let n_reactions = population[best_idx].to_reactions().len();

        println!(
            "gen {:3}  best={:.1}  mean={:.1}  worst={:.1}  reactions={}  genome_len={}",
            generation, best_score, mean_score, worst_score, n_reactions, population[best_idx].len()
        );

        if best_score > best_ever_score {
            best_ever_score = best_score;
            best_ever_genome = Some(Gene8 { genome: population[best_idx].genome.clone(), max_count: population[best_idx].max_count });

            // Save best genome string and reactions JSON
            let genome_str = population[best_idx].to_string();
            let genome_path = cli.output_dir.join(format!("best_gen{:03}.genome", generation));
            std::fs::write(&genome_path, &genome_str).unwrap();

            let reactions_path = cli.output_dir.join("best_reactions.json");
            population[best_idx].save_to_json(reactions_path.to_str().unwrap());
            println!("  → new best! saved to {}", cli.output_dir.display());
        }

        // Build next generation via tournament + crossover + mutation
        // Keep the best individual (elitism)
        let mut next_gen = vec![Gene8 { genome: population[best_idx].genome.clone(), max_count: population[best_idx].max_count }];

        while next_gen.len() < cli.population {
            let p1 = tournament(&population, &scores, 3, &mut rng);
            let p2 = tournament(&population, &scores, 3, &mut rng);
            let (c1, c2) = p1.crossing(p2, &mut rng);
            let c1 = c1.point_mutations(cli.mutation_rate, &mut rng);
            let c2 = c2.point_mutations(cli.mutation_rate, &mut rng);
            next_gen.push(c1);
            if next_gen.len() < cli.population {
                next_gen.push(c2);
            }
        }
        population = next_gen;
    }

    if let Some(best) = &best_ever_genome {
        let final_path = cli.output_dir.join("best_final.json");
        best.save_to_json(final_path.to_str().unwrap());
        println!("\nFinal best chemistry → {}", final_path.display());
        println!("Run it with:");
        let atom_args: Vec<String> = cli.atoms.iter().map(|a| format!("-a {}", a)).collect();
        println!(
            "  cargo run --release --bin simulate -- -r {} --tmax 10000 {}",
            final_path.display(), atom_args.join(" ")
        );
    }
}
