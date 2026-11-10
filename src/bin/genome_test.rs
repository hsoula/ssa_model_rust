use clap::Parser;
use rand::{SeedableRng, rngs::StdRng};
use ssa_model::Gene8::Gene8;
use ssa_model::Reaction;
use std::io::Write;

#[derive(Parser)]
#[command(name = "genome_test", about = "Decode a Gene8 genome and print its reactions")]
struct Cli {
    #[arg(short = 'n', long, default_value = "200")]
    length: usize,

    #[arg(short, long, default_value = "42")]
    seed: u64,

    /// Max value for form counters (states cap at 2×); 10 matches the 10 named AtomTypes
    #[arg(short = 'm', long, default_value = "10")]
    max_count: u8,

    #[arg(short, long, default_value = "data/genome_test")]
    output_dir: std::path::PathBuf,
}

fn main() {
    let cli = Cli::parse();
    let mut rng = StdRng::seed_from_u64(cli.seed);
    let g = Gene8::new_with_max(cli.length, cli.max_count, &mut rng);

    // Pair each pattern with its reaction, keeping only valid ones
    let valid: Vec<(String, String)> = g.extract_patterns_efficient()
        .iter()
        .enumerate()
        .filter_map(|(i, pattern)| {
            let r = Reaction::new_from_vec8(i as u32, pattern.to_vec(), g.max_count);
            if (r.get_substrate_a() >> 8) > 0 && (r.get_substrate_b() >> 8) > 0 {
                let hex = pattern.iter().map(|b| format!("{:X}", b)).collect();
                Some((hex, r.to_string()))
            } else {
                None
            }
        })
        .collect();

    let mut lines: Vec<String> = Vec::new();
    lines.push(format!("genome (seed={} len={}): {}", cli.seed, cli.length, g.to_string()));
    lines.push(String::new());
    lines.push(format!("{} valid reactions:", valid.len()));
    lines.push(String::new());
    for (id, (hex, reaction)) in valid.iter().enumerate() {
        lines.push(format!("[{:3}]  {}  →  {}", id, hex, reaction));
    }

    let output = lines.join("\n");
    println!("{}", output);

    std::fs::create_dir_all(&cli.output_dir).unwrap();
    let filename = cli.output_dir.join(format!("genome_s{}_n{}.txt", cli.seed, cli.length));
    let mut f = std::fs::File::create(&filename).unwrap();
    writeln!(f, "{}", output).unwrap();
    eprintln!("saved → {}", filename.display());
}
