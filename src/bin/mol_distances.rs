use clap::Parser;
use ssa_model::molecules::Molecule;
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "mol_distances",
    about = "Compute pairwise molecule distance matrices from a JSON directory"
)]
struct Cli {
    /// Directory containing molecule JSON files
    #[arg(short, long, default_value = "molecules/json")]
    dir: PathBuf,

    /// Maximum number of molecules to load (0 = all)
    #[arg(long, default_value = "200")]
    max_mols: usize,

    /// WL kernel depth
    #[arg(long, default_value = "2")]
    wl_depth: usize,

    /// Output directory for CSV files
    #[arg(short, long, default_value = "data/distances")]
    output_dir: PathBuf,
}

fn main() {
    let cli = Cli::parse();
    std::fs::create_dir_all(&cli.output_dir).unwrap();

    println!("Loading molecules from {} …", cli.dir.display());
    let mut mols = Molecule::load_from_json_dir(cli.dir.to_str().unwrap());
    println!("  {} loaded", mols.len());

    if cli.max_mols > 0 && mols.len() > cli.max_mols {
        mols.truncate(cli.max_mols);
        println!("  truncated to {}", cli.max_mols);
    }

    let n = mols.len();
    if n == 0 {
        eprintln!("No molecules found – aborting.");
        return;
    }

    // Collect names for CSV header
    let names: Vec<&str> = mols.iter().map(|m| m.name.as_str()).collect();

    // Metrics to compute: (label, closure)
    type DistFn = Box<dyn Fn(&Molecule, &Molecule) -> f64>;
    let wl_depth = cli.wl_depth;
    let metrics: Vec<(&str, DistFn)> = vec![
        ("composition",    Box::new(|a: &Molecule, b: &Molecule| a.composition_distance(b))),
        ("degree_sequence",Box::new(|a: &Molecule, b: &Molecule| a.degree_sequence_distance(b))),
        ("wl",             Box::new(move |a: &Molecule, b: &Molecule| a.wl_distance(b, wl_depth))),
    ];

    for (label, dist_fn) in &metrics {
        println!("Computing {} distances ({n}×{n}) …", label);

        // Build full NxN matrix (symmetric — compute upper triangle, mirror)
        let mut matrix = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let d = dist_fn(&mols[i], &mols[j]);
                matrix[i][j] = d;
                matrix[j][i] = d;
            }
        }

        // Write CSV
        let csv_path = cli.output_dir.join(format!("{}_distance.csv", label));
        let mut f = std::fs::File::create(&csv_path).unwrap();

        // Header row
        write!(f, "name").unwrap();
        for name in &names { write!(f, ",{}", name).unwrap(); }
        writeln!(f).unwrap();

        // Data rows
        for (i, row) in matrix.iter().enumerate() {
            write!(f, "{}", names[i]).unwrap();
            for val in row { write!(f, ",{:.6}", val).unwrap(); }
            writeln!(f).unwrap();
        }

        println!("  → {}", csv_path.display());
    }

    // Summary statistics per metric
    let summary_path = cli.output_dir.join("summary.txt");
    let mut sf = std::fs::File::create(&summary_path).unwrap();
    writeln!(sf, "Molecule distance summary  ({} molecules)", n).unwrap();
    writeln!(sf, "{}", "=".repeat(50)).unwrap();

    for (label, dist_fn) in &metrics {
        let mut vals: Vec<f64> = Vec::with_capacity(n * (n - 1) / 2);
        for i in 0..n {
            for j in (i + 1)..n {
                vals.push(dist_fn(&mols[i], &mols[j]));
            }
        }
        if vals.is_empty() { continue; }
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mean = vals.iter().sum::<f64>() / vals.len() as f64;
        let median = vals[vals.len() / 2];
        writeln!(sf, "\n{label}:  min={:.4}  median={:.4}  mean={:.4}  max={:.4}",
            vals[0], median, mean, vals[vals.len() - 1]).unwrap();
    }
    println!("Summary → {}", summary_path.display());
}
