use clap::Parser;
use std::path::PathBuf;
use ssa_model::reactions::Reaction;

#[derive(Parser)]
#[command(name = "reactions_to_md", about = "Convert a reactions JSON file to readable markdown")]
struct Cli {
    /// Reactions JSON file (produced by Gene8::save_to_json or evolve)
    #[arg(short, long)]
    input: PathBuf,

    /// Write output to file instead of stdout
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();

    let content = std::fs::read_to_string(&cli.input).unwrap_or_else(|e| {
        eprintln!("Cannot read {}: {e}", cli.input.display());
        std::process::exit(1);
    });

    let raw: Vec<serde_json::Value> = serde_json::from_str(&content).unwrap_or_else(|e| {
        eprintln!("Cannot parse JSON: {e}");
        std::process::exit(1);
    });

    let reactions: Vec<Reaction> = raw.iter().enumerate().map(|(i, v)| {
        let parse_key = |field: &str| -> u16 {
            let s = v[field].as_str().unwrap_or("0x0000");
            u16::from_str_radix(s.trim_start_matches("0x").trim_start_matches("0X"), 16)
                .unwrap_or(0)
        };
        let parse_type = |field: &str| -> u8 {
            let s = v[field].as_str().unwrap_or("0b0000");
            u8::from_str_radix(s.trim_start_matches("0b"), 2).unwrap_or(0)
        };
        let s1 = parse_key("s1");
        let s2 = parse_key("s2");
        let p1 = parse_key("p1");
        let p2 = parse_key("p2");
        let sk = ((s1 as u32) << 16) | (s2 as u32);
        let pk = ((p1 as u32) << 16) | (p2 as u32);
        let rt = parse_type("reaction_type");
        let rs = v["reaction_speed"].as_f64().unwrap_or(0.1);
        Reaction::new(i as u32)
            .set_substrate_key(sk)
            .set_product_key(pk)
            .set_reaction_type(rt)
            .set_reaction_speed(rs)
    }).collect();

    let mut md = String::new();
    md.push_str(&format!("# Reactions  ({})\n\n", cli.input.display()));
    md.push_str(&format!("Loaded {} reactions.\n\n", reactions.len()));
    md.push_str("| # | Reaction | Type | Speed |\n");
    md.push_str("|---|----------|------|-------|\n");

    let type_label = |rt: u8| match rt & 0b11 {
        0b00 => "contact",
        0b01 => "→ bond",
        0b10 => "break bond",
        0b11 => "exchange",
        _    => "?",
    };

    for (i, r) in reactions.iter().enumerate() {
        md.push_str(&format!("| {:>3} | `{}` | {} | {:.2} |\n",
            i + 1,
            r.to_string(),
            type_label(r.reaction_type()),
            r.reaction_speed(),
        ));
    }

    match cli.output {
        Some(path) => {
            std::fs::write(&path, &md).unwrap_or_else(|e| {
                eprintln!("Cannot write {}: {e}", path.display());
                std::process::exit(1);
            });
            eprintln!("Written to {}", path.display());
        }
        None => print!("{}", md),
    }
}
