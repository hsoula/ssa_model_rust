/// gene_doc — print the Gene8 encoding specification + worked examples.
///
/// Run this after any change to the token table in Gene8.rs / reactions.rs to
/// regenerate docs/gene8_encoding.md:
///
///   cargo run --release --bin gene_doc > docs/gene8_encoding.md

use rand::{SeedableRng, rngs::StdRng};
use ssa_model::Gene8::Gene8;

// ── Token table (single source of truth for the doc) ─────────────────────────

struct Token {
    hex: &'static str,
    role: &'static str,
    field: &'static str,
}

const TOKENS: &[Token] = &[
    Token { hex: "E (0x0E)", role: "START marker", field: "opens a reaction pattern"              },
    Token { hex: "F (0x0F)", role: "END marker",   field: "closes a reaction pattern"             },
    Token { hex: "0,1",      role: "sf1++",        field: "substrate 1 form  (x2 probability)"   },
    Token { hex: "2,3",      role: "ss1++",        field: "substrate 1 state (x2 probability)"   },
    Token { hex: "4,5",      role: "sf2++",        field: "substrate 2 form  (x2 probability)"   },
    Token { hex: "6,7",      role: "ss2++",        field: "substrate 2 state (x2 probability)"   },
    Token { hex: "8",        role: "contact++",    field: "reaction-type bits (value mod 4)"      },
    Token { hex: "9,A",      role: "ps1++",        field: "product 1 state   (x2 probability)"   },
    Token { hex: "B,C",      role: "ps2++",        field: "product 2 state   (x2 probability)"   },
    Token { hex: "D",        role: "speed += 0.1", field: "reaction speed increment"              },
];

// ── Helpers ───────────────────────────────────────────────────────────────────

fn reaction_type_name(t: u8) -> &'static str {
    match t & 0b11 {
        0b00 => "contact (no bond change)",
        0b01 => "product bond",
        0b10 => "substrate bond",
        0b11 => "bond exchange",
        _    => "unknown",
    }
}

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    println!("# Gene8 Encoding");
    println!();
    println!("Gene8 is a linear genome of u8 tokens, each in `[0x00, 0x0F]` (one hex digit).");
    println!("Reactions are encoded as sub-sequences delimited by **E** (start) and **F** (end).");
    println!("The genome is created by sampling uniformly from 16 symbols and reducing mod 16.");
    println!();

    // ── Token table ──────────────────────────────────────────────────────────
    println!("## Token table");
    println!();
    println!("| Token | Role | Field |");
    println!("|-------|------|-------|");
    for t in TOKENS {
        println!("| `{}` | {} | {} |", t.hex, t.role, t.field);
    }
    println!();

    // ── Decoding rules ───────────────────────────────────────────────────────
    println!("## Decoding rules");
    println!();
    println!("For each `E…F` pattern the payload tokens are tallied into counters:");
    println!();
    println!("```");
    println!("sf1, ss1, sf2, ss2   ← substrate forms/states");
    println!("ps1, ps2             ← product states  (token 5 → ps1, token 6 → ps2)");
    println!("contact              ← reaction-type bits (mod 4)");
    println!("speed                ← 0.1 × (count of token 9)");
    println!("```");
    println!();
    println!("Atom keys are 16-bit words: `(form << 8) | state`.");
    println!();
    println!("```");
    println!("s1 = (sf1 << 8) | ss1      s2 = (sf2 << 8) | ss2");
    println!("p1 = (sf1 << 8) | ps1      p2 = (sf2 << 8) | ps2");
    println!("```");
    println!();
    println!("Product forms **always** equal substrate forms — only states can change.");
    println!("A reaction is **invalid** (discarded) if `sf1 == 0` or `sf2 == 0`.");
    println!();

    // ── Reaction type bits ───────────────────────────────────────────────────
    println!("## Reaction-type bits");
    println!();
    println!("| contact mod 4 | Meaning |");
    println!("|---------------|---------|");
    for v in 0u8..4 {
        println!("| `0b{:02b}` ({}) | {} |", v, v, reaction_type_name(v));
    }
    println!();

    // ── Worked examples ──────────────────────────────────────────────────────
    println!("## Worked examples");
    println!();

    let examples: &[(&str, &[u8])] = &[
        (
            "Minimal bond-forming reaction",
            &[0x0E, 0x00, 0x02, 0x04, 0x0F],
            // E 0 2 4 F
            // sf1=1 sf2=1 contact=1 → s1=(1<<8)|0=0x0100  s2=(1<<8)|0=0x0200
            //                          p1=s1  p2=s2  type=1 (product bond)
        ),
        (
            "State-change reaction with speed",
            &[0x0E, 0x00, 0x01, 0x02, 0x05, 0x09, 0x0F],
            // sf1=1 ss1=1 sf2=1 ps1=1 speed=0.1
            // s1=(1<<8)|1=0x0101  s2=(1<<8)|0=0x0200
            // p1=(1<<8)|1=0x0101  p2=(1<<8)|0=0x0200  speed=0.1
        ),
        (
            "Two-substrate bond-exchange with product state changes",
            &[0x0E, 0x00, 0x02, 0x04, 0x04, 0x04, 0x05, 0x07, 0x09, 0x09, 0x0F],
            // sf1=1 sf2=1 contact=3 ps1=1 ps2=1 speed=0.2
        ),
    ];

    for (title, pattern) in examples {
        let hex: String = pattern.iter().map(|b| format!("{:X}", b)).collect();
        println!("### {}", title);
        println!();
        println!("Pattern: `{}`", hex);
        println!();

        // Tally counters manually to annotate
        let (mut sf1, mut ss1, mut sf2, mut ss2) = (0u8, 0u8, 0u8, 0u8);
        let (mut ps1, mut ps2, mut contact) = (0u8, 0u8, 0u8);
        let mut speed = 0.0f64;
        for &b in *pattern {
            match b {
                0x0E | 0x0F => {}
                0x00 | 0x01 => sf1 += 1,
                0x02 | 0x03 => ss1 += 1,
                0x04 | 0x05 => sf2 += 1,
                0x06 | 0x07 => ss2 += 1,
                0x08        => contact += 1,
                0x09 | 0x0A => ps1 += 1,
                0x0B | 0x0C => ps2 += 1,
                0x0D        => speed += 0.1,
                _ => {}
            }
        }
        contact %= 4;

        println!("| Counter | Value |");
        println!("|---------|-------|");
        println!("| sf1 | {} |", sf1);
        println!("| ss1 | {} |", ss1);
        println!("| sf2 | {} |", sf2);
        println!("| ss2 | {} |", ss2);
        println!("| contact (mod 4) | {} → {} |", contact, reaction_type_name(contact));
        println!("| ps1 | {} |", ps1);
        println!("| ps2 | {} |", ps2);
        println!("| speed | {:.1} |", speed);
        println!();

        // Only show reaction string if valid
        if sf1 > 0 && sf2 > 0 {
            let s1 = ((sf1 as u16) << 8) | (ss1 as u16);
            let s2 = ((sf2 as u16) << 8) | (ss2 as u16);
            let p1 = ((sf1 as u16) << 8) | (ps1 as u16);
            let p2 = ((sf2 as u16) << 8) | (ps2 as u16);

            let sym = |k: u16| -> String {
                let form  = ((k >> 8) as u8 + 96) as char;
                let state = (k & 0xFF) as u8;
                format!("{}{}", form, state)
            };
            let sep_sub = if contact & 0b10 != 0 { "+" } else { "." };
            let sep_prd = if contact & 0b01 != 0 { "+" } else { "." };
            println!("**Reaction:** `{}{}{} → {}{}{}`  *(speed={:.1})*",
                sym(s1), sep_sub, sym(s2), sym(p1), sep_prd, sym(p2), speed);
        } else {
            println!("**Invalid** — sf1={} sf2={}: reaction discarded.", sf1, sf2);
        }
        println!();
    }

    // ── Random genome example ─────────────────────────────────────────────────
    println!("## Random genome example (seed=42, len=60)");
    println!();
    let mut rng = StdRng::seed_from_u64(42);
    let g = Gene8::new(60, &mut rng);
    println!("Genome: `{}`", g.to_string());
    println!();
    let reactions = g.to_reactions();
    println!("{} valid reaction(s) decoded:", reactions.len());
    println!();
    for (i, r) in reactions.iter().enumerate() {
        println!("{}. `{}`", i + 1, r.to_string());
    }
    println!();

    eprintln!("Pipe this output to docs/gene8_encoding.md to update the documentation.");
}
