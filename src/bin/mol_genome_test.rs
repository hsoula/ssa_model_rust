//! mol_genome_test — decode a MolGenome and print the resulting molecules.
//!
//! Usage:
//!   mol_genome_test --genome E100201000F --max-form 2 --max-state 4
//!   mol_genome_test --random --size 120 --max-form 3 --max-state 4 --seed 42
//!   mol_genome_test --random --size 200 --max-form 2 --max-state 2 --save out.json

use clap::Parser;
use rand::{SeedableRng, rngs::StdRng};
use ssa_model::mol_genome::MolGenome;

#[derive(Parser)]
#[command(name = "mol_genome_test",
          about = "Decode a MolGenome and display the molecules it encodes")]
struct Cli {
    /// Hex genome string (e.g. "E100201000F").  Mutually exclusive with --random.
    #[arg(long)]
    genome: Option<String>,

    /// Generate a random genome instead of using --genome
    #[arg(long, default_value = "false")]
    random: bool,

    /// Length of the random genome (tokens)
    #[arg(long, default_value = "120")]
    size: usize,

    /// Maximum atom form (forms are 1..=max_form)
    #[arg(long, default_value = "2")]
    max_form: u8,

    /// Maximum atom state (states are 0..=max_state)
    #[arg(long, default_value = "4")]
    max_state: u8,

    /// Random seed
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Save decoded genome+molecules to this JSON file
    #[arg(long)]
    save: Option<String>,
}

fn atom_sym(form: u8, state: u8) -> String {
    format!("{}{}", (form + 96) as char, state)
}

fn main() {
    let cli = Cli::parse();
    let mut rng = StdRng::seed_from_u64(cli.seed);

    // ── Build genome ──────────────────────────────────────────────────────────
    let gene = if cli.random || cli.genome.is_none() {
        println!("Generating random genome  size={}  max_form={}  max_state={}  seed={}",
            cli.size, cli.max_form, cli.max_state, cli.seed);
        MolGenome::new(cli.size, cli.max_form, cli.max_state, &mut rng)
    } else {
        let s = cli.genome.unwrap();
        println!("Decoding genome: {s}");
        MolGenome::from_hex(&s, cli.max_form, cli.max_state)
    };

    // ── Genome summary ────────────────────────────────────────────────────────
    println!("\n══ Genome ══════════════════════════════════════════════════════");
    println!("  Hex : {}", gene.to_hex());
    println!("  Len : {} tokens", gene.len());

    let n_start = gene.genome.iter().filter(|&&b| b == 0x0E).count();
    let n_end   = gene.genome.iter().filter(|&&b| b == 0x0F).count();
    let n_state = gene.genome.iter().filter(|&&b| b == 0x00).count();
    let n_form  = gene.genome.iter().filter(|&&b| b > 0x00 && b < 0x0E).count();
    println!("  E (start) : {n_start}  F (end) : {n_end}  \
              0 (state) : {n_state}  form tokens : {n_form}");

    // ── Patterns ──────────────────────────────────────────────────────────────
    let patterns = gene.extract_patterns();
    println!("\n══ Patterns ({} found) ══════════════════════════════════════════", patterns.len());
    for (i, pat) in patterns.iter().enumerate() {
        let pat_str: String = pat.iter().map(|b| format!("{:X}", b)).collect();
        println!("  [{i}] {pat_str}  ({} tokens)", pat.len());
    }

    // ── Molecules ─────────────────────────────────────────────────────────────
    let mols = gene.to_molecules();
    println!("\n══ Decoded molecules ({} valid, ≥2 atoms) ═══════════════════════", mols.len());

    if mols.is_empty() {
        println!("  (none — no pattern contained ≥ 2 atoms)");
    }

    for (i, m) in mols.iter().enumerate() {
        println!("\n  ── Molecule [{}]  {} atoms, {} bonds ──────────────────────────",
            i, m.atoms.len(), m.links.len());

        // Atoms sorted by id
        let mut atoms: Vec<(u32, u8, u8)> = m.atoms.iter()
            .map(|(&id, &(f, s))| (id, f, s))
            .collect();
        atoms.sort_by_key(|(id, _, _)| *id);

        let chain: Vec<String> = atoms.iter()
            .map(|(_, f, s)| atom_sym(*f, *s))
            .collect();
        println!("    Chain  : {}", chain.join(" — "));

        println!("    Atoms  :");
        for (id, f, s) in &atoms {
            println!("      id={id}  form={}  state={s}  ({})",
                f, atom_sym(*f, *s));
        }

        println!("    Bonds  :");
        let mut bonds: Vec<(u32, u32, u32)> = m.links.iter()
            .map(|(&id, &(a, b))| (id, a, b))
            .collect();
        bonds.sort_by_key(|(id, _, _)| *id);
        for (id, a, b) in &bonds {
            let sa = atoms.iter().find(|(i,_,_)| i==a).map(|(_,f,s)| atom_sym(*f,*s)).unwrap_or_default();
            let sb = atoms.iter().find(|(i,_,_)| i==b).map(|(_,f,s)| atom_sym(*f,*s)).unwrap_or_default();
            println!("      bond {id}: {a}({sa}) — {b}({sb})");
        }
    }

    // ── Composition summary across all molecules ──────────────────────────────
    if !mols.is_empty() {
        use std::collections::HashMap;
        let mut comp: HashMap<(u8,u8), usize> = HashMap::new();
        for m in &mols {
            for &(f, s) in m.atoms.values() {
                *comp.entry((f, s)).or_insert(0) += 1;
            }
        }
        let mut comp_v: Vec<_> = comp.iter().collect();
        comp_v.sort_by_key(|((f,s),_)| (*f, *s));
        println!("\n══ Atom composition across all molecules ════════════════════════");
        for ((f, s), n) in &comp_v {
            println!("  {}  ×{}", atom_sym(*f, *s), n);
        }
        let total_atoms: usize = mols.iter().map(|m| m.atoms.len()).sum();
        let total_bonds: usize = mols.iter().map(|m| m.links.len()).sum();
        println!("  Total: {total_atoms} atoms, {total_bonds} bonds");
    }

    // ── Optional JSON save ────────────────────────────────────────────────────
    if let Some(path) = cli.save {
        gene.save_to_json(&path);
        println!("\nSaved → {path}");
    }
}
