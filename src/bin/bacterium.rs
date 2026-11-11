//! Demo: one bacterium in a medium fed with A1 (uppercase = outside, lowercase = inside).
//!
//!   cargo run --release --bin bacterium -- --tmax 200 --dt 20

use clap::Parser;
use rand::{SeedableRng, rngs::StdRng};
use ssa_model::bacterium::{Bacterium, BacteriumReactor, Gene};

#[derive(Parser)]
#[command(name = "bacterium", about = "Gillespie run of a single bacterium in a medium")]
struct Cli {
    #[arg(long, default_value = "42")]
    seed: u64,
    #[arg(long, default_value = "200")]
    tmax: f64,
    /// report interval
    #[arg(long, default_value = "20")]
    dt: f64,
}

fn genome() -> Result<Bacterium, String> {
    let mut b = Bacterium::new();
    // passive channel: A1 in, a1 out (dragging whatever molecule it is bonded to)
    b.add_gene(Gene::new("permease").with("A1 -> a1", 0.005)?);
    b.add_gene(Gene::new("efflux").with("a1 -> A1", 0.01)?);
    // receptor b2 grabs an outside A1, pulls it in (a2) and becomes active (b3), then releases
    b.add_gene(
        Gene::new("receptor")
            .with("A1 + b2 -> a2.b3", 0.002)?
            .with("a2.b3 -> a2 + b2", 0.5)?,
    );
    Ok(b)
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    let mut rng = StdRng::seed_from_u64(cli.seed);

    let mut br = BacteriumReactor::new(genome()?);
    br.add_environment_rule("0 -> A1", 1.0)?; // food supply
    br.add_environment_rule("A1 -> 0", 0.01)?; // dilution
    br.add_atoms("A1", 50)?;
    br.add_chains("C4.A1", 10)?; // A1 carrying a C4: C4 comes in with it
    br.add_atoms("b2", 20)?;

    print!("{}", br.report());
    let mut t = 0.0;
    while t < cli.tmax {
        t += cli.dt;
        br.run(t, &mut rng);
        print!("{}", br.report());
    }
    Ok(())
}
