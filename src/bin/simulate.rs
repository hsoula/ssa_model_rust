use clap::Parser;
use rand::{SeedableRng, rngs::StdRng};
use ssa_model::reactions::Reaction;
use ssa_model::reactor::Reactor;
use ssa_model::molecules::Molecule;
use ssa_model::utils::parse_atom_spec;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "simulate", about = "Run an SSA chemistry simulation")]
struct Cli {
    /// Reaction JSON files (repeat for multiple: -r a.json -r b.json)
    #[arg(short, long = "reactions", required = true, num_args = 1..)]
    reactions: Vec<PathBuf>,

    /// Initial atoms as KEY:COUNT (hex key, e.g. 0x0100:300); repeat for multiple types.
    /// Either byte of the key may contain a '*' to wildcard it, e.g. 0x0*0*:50 (every
    /// form × every state), 0x0*03:50 (every form, state=3), 0x01**:50 (form=1, every state).
    /// May be omitted entirely if --mol supplies the whole initial population.
    #[arg(short, long = "atoms", num_args = 1..)]
    atoms: Vec<String>,

    /// Max form used to expand a wildcarded (*) form in an atom spec (forms 1..=max_form)
    #[arg(long, default_value = "10")]
    max_form: u8,

    /// Max state used to expand a wildcarded (*) state in an atom spec (states 0..=max_state)
    #[arg(long, default_value = "10")]
    max_state: u8,

    /// Random seed
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Max simulation time
    #[arg(long, default_value = "10000")]
    tmax: f64,

    /// Output directory for molecule files and summary
    #[arg(short, long = "output-dir", default_value = "data/output")]
    output_dir: PathBuf,

    /// Minimum molecule size to save
    #[arg(long, default_value = "3")]
    min_size: usize,

    /// Load an initial molecule from a .mol file
    #[arg(long)]
    mol: Option<PathBuf>,

    /// Skip writing output files
    #[arg(long)]
    no_save: bool,
}

fn main() {
    let cli = Cli::parse();
    if cli.atoms.is_empty() && cli.mol.is_none() {
        eprintln!("Error: provide at least one of -a KEY:COUNT or --mol FILE");
        std::process::exit(1);
    }
    let mut rng = StdRng::seed_from_u64(cli.seed);
    let mut reactor = Reactor::new();

    // Load reactions
    for path in &cli.reactions {
        let reactions = Reaction::news_from_json(path.to_str().unwrap(), 1);
        println!("Loaded {} reactions from {}", reactions.len(), path.display());
        for r in reactions {
            reactor.add_reaction(r);
        }
    }

    // Load optional initial molecule
    if let Some(mol_path) = &cli.mol {
        let mut molecule = Molecule::new();
        molecule.load_from_file(mol_path.to_str().unwrap());
        reactor.add_molecule(&molecule);
        println!("Loaded molecule from {}", mol_path.display());
    }

    // Add initial atoms
    for spec in &cli.atoms {
        let pairs = parse_atom_spec(spec, cli.max_form, cli.max_state).unwrap_or_else(|e| {
            eprintln!("Error parsing atom spec '{spec}': {e}");
            std::process::exit(1);
        });
        for (key, count) in pairs {
            reactor.add_atom_from_keys(key, count);
            println!("Added {} atoms of key 0x{:04X}", count, key);
        }
    }

    println!("\n=== Initial state ===");
    reactor.display_stats();

    let t_final = reactor.main_loop(cli.tmax, &mut rng);

    println!("\n=== Final state  t = {:.3} ===", t_final);
    reactor.display_events();
    let (atom_stats, link_stats) = reactor.display_stats();

    let molecules = reactor.find_molecules();
    let size_dist: HashMap<usize, usize> = molecules.iter()
        .map(|m| m.size())
        .fold(HashMap::new(), |mut acc, s| { *acc.entry(s).or_insert(0) += 1; acc });

    println!("\n=== Molecule size distribution ===");
    let mut sizes: Vec<_> = size_dist.iter().collect();
    sizes.sort_by_key(|(k, _)| *k);
    for (size, count) in &sizes {
        println!("  size {:4}: {:4} molecule(s)", size, count);
    }

    if !cli.no_save {
        std::fs::create_dir_all(&cli.output_dir).unwrap();

        for (i, mol) in molecules.iter().enumerate() {
            if mol.size() >= cli.min_size {
                let path = cli.output_dir.join(format!("mol_{:04}.txt", i));
                mol.save_to_file(path.to_str().unwrap());
            }
        }

        let summary = serde_json::json!({
            "seed": cli.seed,
            "tmax": cli.tmax,
            "t_final": t_final,
            "reaction_files": cli.reactions.iter().map(|p| p.to_string_lossy()).collect::<Vec<_>>(),
            "initial_atoms": cli.atoms,
            "n_molecules_total": molecules.len(),
            "size_distribution": size_dist.iter()
                .map(|(k, v)| (k.to_string(), *v))
                .collect::<HashMap<_, _>>(),
            "atom_stats": atom_stats,
            "link_stats": link_stats.iter()
                .map(|(k, v)| (format!("0x{:08X}", k), *v))
                .collect::<HashMap<_, _>>(),
        });

        let summary_path = cli.output_dir.join("summary.json");
        std::fs::write(&summary_path, serde_json::to_string_pretty(&summary).unwrap()).unwrap();
        println!("\nOutput written to {}", cli.output_dir.display());
    }
}
