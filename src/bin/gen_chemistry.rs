use clap::Parser;
use rand::{SeedableRng, rngs::StdRng, Rng};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "gen_chemistry", about = "Generate a random reaction set as JSON")]
struct Cli {
    /// Number of reactions to generate
    #[arg(short = 'n', long, default_value = "10")]
    n_reactions: usize,

    /// Number of distinct atom forms (types), 1..n_forms
    #[arg(short = 'f', long, default_value = "4")]
    n_forms: u8,

    /// Number of distinct atom states per form, 0..n_states
    #[arg(short = 's', long, default_value = "4")]
    n_states: u8,

    /// Include production reactions (0 → atom) at ~10% probability
    #[arg(long,default_value = "false")]
    production: bool,

    /// Include degradation reactions (atom → 0) at ~10% probability
    #[arg(long,default_value = "false")]
    degradation: bool,

    /// Fraction of reactions that bond substrates (substrate-contact type)
    #[arg(long, default_value = "0.5")]
    contact_frac: f64,

    /// Min reaction speed
    #[arg(long, default_value = "0.01")]
    speed_min: f64,

    /// Max reaction speed
    #[arg(long, default_value = "1.0")]
    speed_max: f64,

    /// Random seed
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Output JSON file
    #[arg(short, long, default_value = "data/chemistry.json")]
    output: PathBuf,
}

// ── Helpers for markdown output ───────────────────────────────────────────────
use ssa_model::reactions::Reaction;

/// Convert an atom u16 key (form<<8 | state) to its display symbol.
/// Production source / degradation sink are represented as "0".
fn atom_sym(k: u16) -> String {
    if k == 0 { return "0".to_string(); }
    let form  = ((k >> 8) & 0xff) as u8;
    let state = (k & 0xff) as u8;
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

/// Sort key: (bucket, sub_a, sub_b)
/// bucket 0 = production (sub is "0"), 1 = normal, 2 = degradation (prod is "0")
fn sort_key(r: &Reaction) -> (u8, u16, u16) {
    if r.is_production()  { (0, r.get_product_a(),  0) }
    else if r.is_degradation() { (2, r.get_substrate_a(), 0) }
    else { (1, r.get_substrate_a(), r.get_substrate_b()) }
}

fn generate_markdown(rs: &[Reaction], json_path: &std::path::Path) -> String {
    use std::fmt::Write;
    let mut md = String::new();
    let _ = writeln!(md, "# Chemistry — `{}`\n", json_path.display());
    let _ = writeln!(md, "**{} reactions**\n", rs.len());

    // ── Sorted reaction table ─────────────────────────────────────────────────
    let _ = writeln!(md, "| # | Sub A | Sub B | Prod A | Prod B | Type | Speed |");
    let _ = writeln!(md, "|---|-------|-------|--------|--------|------|-------|");

    let mut sorted: Vec<&Reaction> = rs.iter().collect();
    sorted.sort_by_key(|r| sort_key(r));

    for (i, r) in sorted.iter().enumerate() {
        let (sub_a, sub_b, prod_a, prod_b) = if r.is_production() {
            ("0".to_string(), "—".to_string(),
             atom_sym(r.get_product_a()), "—".to_string())
        } else if r.is_degradation() {
            (atom_sym(r.get_substrate_a()), "—".to_string(),
             "0".to_string(), "—".to_string())
        } else {
            (atom_sym(r.get_substrate_a()), atom_sym(r.get_substrate_b()),
             atom_sym(r.get_product_a()),   atom_sym(r.get_product_b()))
        };
        let _ = writeln!(md, "| {} | `{}` | `{}` | `{}` | `{}` | {} | {:.3} |",
            i + 1, sub_a, sub_b, prod_a, prod_b,
            rxn_type_label(r), r.reaction_speed());
    }
    let _ = writeln!(md);

    // ── Per-substrate outcome matrix (who reacts with whom → what) ────────────
    // Collect all unique atom keys that appear as substrate A
    let normal: Vec<&&Reaction> = sorted.iter()
        .filter(|r| !r.is_production() && !r.is_degradation())
        .collect();

    if !normal.is_empty() {
        // All distinct substrate A values, sorted by (form, state)
        let mut sub_a_keys: Vec<u16> = normal.iter()
            .map(|r| r.get_substrate_a()).collect();
        sub_a_keys.sort();
        sub_a_keys.dedup();

        let _ = writeln!(md, "## Substrate A outcomes\n");
        let _ = writeln!(md, "Rows = Substrate A, sorted by form then state.\n");
        let _ = writeln!(md, "| Sub A | Sub B | → | Prod A | Prod B | Type | Speed |");
        let _ = writeln!(md, "|-------|-------|---|--------|--------|------|-------|");

        for sa in &sub_a_keys {
            let rows: Vec<&&&Reaction> = normal.iter()
                .filter(|r| r.get_substrate_a() == *sa)
                .collect();
            for (j, r) in rows.iter().enumerate() {
                let label = format!("`{}`", atom_sym(*sa));
                let _ = writeln!(md, "| {} | `{}` | → | `{}` | `{}` | {} | {:.3} |",
                    label,
                    atom_sym(r.get_substrate_b()),
                    atom_sym(r.get_product_a()),
                    atom_sym(r.get_product_b()),
                    rxn_type_label(r), r.reaction_speed());
            }
        }
        let _ = writeln!(md);
    }

    md
}

fn main() {
    let cli = Cli::parse();
    let mut rng = StdRng::seed_from_u64(cli.seed);

    let mut reactions: Vec<serde_json::Value> = Vec::new();
    let reac_types = vec!["0b0000", "0b0001", "0b0010", "0b0011"];
    for _ in 0..cli.n_reactions {
        let reaction_type = if cli.production && rng.random_bool(0.1) {
            "0b0100"  // production
        } else if cli.degradation && rng.random_bool(0.1) {
            "0b1000"  // degradation
        } else {
            reac_types[rng.random_range(0..4)]
        };

        let sf1: u16 = rng.random_range(1u16..=cli.n_forms as u16);
        let ss1: u16 = rng.random_range(0u16..cli.n_states as u16);
        let sf2: u16 = rng.random_range(1u16..=cli.n_forms as u16);
        let ss2: u16 = rng.random_range(0u16..cli.n_states as u16);
        let pf1: u16 = sf1;
        let ps1: u16 = rng.random_range(0u16..cli.n_states as u16);
        let pf2: u16 = sf2;
        let ps2: u16 = rng.random_range(0u16..cli.n_states as u16);
        let speed: f64 = rng.random_range(cli.speed_min..cli.speed_max);

        let s1 = (sf1 << 8) | ss1;
        let s2 = (sf2 << 8) | ss2;
        let p1 = (pf1 << 8) | ps1;
        let p2 = (pf2 << 8) | ps2;

        reactions.push(serde_json::json!({
            "s1": format!("0x{:04X}", s1),
            "s2": format!("0x{:04X}", s2),
            "p1": format!("0x{:04X}", p1),
            "p2": format!("0x{:04X}", p2),
            "reaction_type": reaction_type,
            "reaction_speed": (speed * 1000.0).round() / 1000.0,
        }));
    }

    if let Some(parent) = cli.output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).unwrap();
        }
    }
    std::fs::write(&cli.output, serde_json::to_string_pretty(&reactions).unwrap()).unwrap();
    println!("Generated {} reactions → {}", reactions.len(), cli.output.display());

    // Load back through Reaction for human-readable output and markdown
    let rs = Reaction::news_from_json(cli.output.to_str().unwrap(), 1);
    for r in &rs {
        println!("  {}", r.to_string());
    }

    // Write markdown report alongside the JSON (same stem, .md extension)
    let md_path = cli.output.with_extension("md");
    let md = generate_markdown(&rs, &cli.output);
    std::fs::write(&md_path, &md).unwrap_or_else(|e| {
        eprintln!("Cannot write markdown {}: {e}", md_path.display());
    });
    println!("Markdown report → {}", md_path.display());
}
