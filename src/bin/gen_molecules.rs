//! gen_molecules — generate a batch of random molecules and save them to JSON + mol format.
//!
//! Each molecule is produced by `Molecule::random()`:
//!   • Connected via a random spanning tree (Prüfer sequence).
//!   • Atom forms sampled uniformly from 1..=n_forms; states always start at 0.
//!   • Optional extra edges add cycles.
//!
//! Outputs (same base path, different extensions):
//!   <output>.json  — JSON array: [{name, n_atoms, n_bonds, atoms:[{id,form,state}…], bonds:[{id,a,b}…]}…]
//!   <output>.mol   — concatenated mol blocks, one per molecule, each terminated by "//"
//!
//! The .mol format per molecule:
//!   ATOM <n>
//!   <id> <form> <state>
//!   …
//!   BONDS
//!   <atom_a> <atom_b>
//!   …
//!   //
//!
//! Usage examples:
//!   gen_molecules -n 50 -o data/mols
//!   gen_molecules -n 20 --min-atoms 3 --max-atoms 8 -f 3 --extra-edges 1 --seed 7 -o data/mols

use clap::Parser;
use rand::{SeedableRng, rngs::StdRng, Rng};
use ssa_model::molecules::Molecule;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "gen_molecules",
          about = "Generate random connected molecules and save as JSON + mol")]
struct Cli {
    /// Number of molecules to generate
    #[arg(short = 'n', long, default_value = "20")]
    n_molecules: usize,

    /// Minimum number of atoms per molecule
    #[arg(long, default_value = "2")]
    min_atoms: usize,

    /// Maximum number of atoms per molecule
    #[arg(long, default_value = "6")]
    max_atoms: usize,

    /// Number of distinct atom forms (types), forms are 1..=n_forms
    #[arg(short = 'f', long, default_value = "3")]
    n_forms: u8,

    /// Number of distinct atom forms (types), forms are 1..=n_forms
    #[arg(short = 's', long, default_value = "3")]
    n_states: u8,
    /// Number of extra edges added beyond the spanning tree (adds cycles)
    #[arg(long, default_value = "0")]
    extra_edges: usize,

    /// Random seed
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Output base path (extensions .json and .mol are appended automatically)
    #[arg(short, long, default_value = "data/molecules")]
    output: String,
}

// ── Serialisers ───────────────────────────────────────────────────────────────

fn mol_to_json(mol: &Molecule) -> serde_json::Value {
    let mut atoms: Vec<serde_json::Value> = mol.atoms.iter()
        .map(|(&id, &(form, state))| serde_json::json!({
            "id": id, "form": form, "state": state
        }))
        .collect();
    atoms.sort_by_key(|v| v["id"].as_u64().unwrap_or(0));

    let mut bonds: Vec<serde_json::Value> = mol.links.iter()
        .map(|(&id, &(a, b))| serde_json::json!({
            "id": id, "a": a, "b": b
        }))
        .collect();
    bonds.sort_by_key(|v| v["id"].as_u64().unwrap_or(0));

    serde_json::json!({
        "name":    mol.name,
        "n_atoms": mol.atoms.len(),
        "n_bonds": mol.links.len(),
        "atoms":   atoms,
        "bonds":   bonds,
    })
}

/// Render one molecule as a mol-format block (ATOM … BONDS … //).
fn mol_to_mol_block(mol: &Molecule) -> String {
    let mut out = String::new();

    // Optional name comment so the file is self-documenting
    out.push_str(&format!("# {}\n", mol.name));

    out.push_str(&format!("ATOM {}\n", mol.atoms.len()));
    let mut atom_ids: Vec<u32> = mol.atoms.keys().cloned().collect();
    atom_ids.sort();
    for id in &atom_ids {
        let (form, state) = mol.atoms[id];
        out.push_str(&format!("{id} {form} {state}\n"));
    }

    out.push_str("BONDS\n");
    let mut link_ids: Vec<u32> = mol.links.keys().cloned().collect();
    link_ids.sort();
    for lid in &link_ids {
        let (a, b) = mol.links[lid];
        out.push_str(&format!("{a} {b}\n"));
    }

    out.push_str("//\n");
    out
}

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let cli = Cli::parse();
    let mut rng = StdRng::seed_from_u64(cli.seed);

    if cli.min_atoms == 0 || cli.min_atoms > cli.max_atoms {
        eprintln!("Error: need 1 ≤ min_atoms ≤ max_atoms");
        std::process::exit(1);
    }
    if cli.n_forms == 0 {
        eprintln!("Error: n_forms must be ≥ 1");
        std::process::exit(1);
    }

    // Ensure output directory exists (strip any extension from the base path)
    let base = PathBuf::from(&cli.output);
    if let Some(parent) = base.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).unwrap();
        }
    }

    let json_path = base.with_extension("json");
    let mol_path  = base.with_extension("mol");

    let forms: Vec<u8> = (1..=cli.n_forms).collect();
    let states: Vec<u8> = (0..=cli.n_states).collect();

    let mut json_vals: Vec<serde_json::Value> = Vec::with_capacity(cli.n_molecules);
    let mut mol_blocks = String::new();

    for i in 0..cli.n_molecules {
        let n = if cli.min_atoms == cli.max_atoms {
            cli.min_atoms
        } else {
            rng.random_range(cli.min_atoms..=cli.max_atoms)
        };
        let mut mol = Molecule::random(n, &forms, &states, cli.extra_edges, &mut rng);
        mol.name = format!("mol_{i:04}");

        json_vals.push(mol_to_json(&mol));
        mol_blocks.push_str(&mol_to_mol_block(&mol));
    }

    // Write JSON
    std::fs::write(&json_path, serde_json::to_string_pretty(&json_vals).unwrap())
        .unwrap_or_else(|e| { eprintln!("Cannot write {}: {e}", json_path.display()); std::process::exit(1); });

    // Write mol
    std::fs::write(&mol_path, &mol_blocks)
        .unwrap_or_else(|e| { eprintln!("Cannot write {}: {e}", mol_path.display()); std::process::exit(1); });

    // Summary
    let sizes: Vec<usize> = json_vals.iter()
        .map(|v| v["n_atoms"].as_u64().unwrap_or(0) as usize)
        .collect();
    let min_s  = sizes.iter().min().copied().unwrap_or(0);
    let max_s  = sizes.iter().max().copied().unwrap_or(0);
    let mean_s = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;

    println!("Generated {} molecules", json_vals.len());
    println!("  atoms : min={min_s}  max={max_s}  mean={mean_s:.1}  forms=1..={}", cli.n_forms);
    println!("  JSON  → {}", json_path.display());
    println!("  mol   → {}", mol_path.display());
}
