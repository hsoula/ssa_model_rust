//! Bacterium reactor
//!
//! A single bacterium sits in a medium. The compartment of an atom is encoded in the case
//! of its letter (one bit of the form, `atoms::OUTSIDE`):
//!
//! * `A1`, `B3` ...  uppercase  -> outside (medium)
//! * `a1`, `b3` ...  lowercase  -> inside the bacterium
//!
//! Crossing the membrane only flips the case; the state is kept. Molecules cross as a whole:
//! when a reaction moves one atom, every atom bonded to it (directly or not) follows, so with
//! `B1.A1` outside and the enzyme `A1 -> a1`, the result is `b1.a1` inside.
//! As a consequence a molecule always lives on one side, and a rule may not bond two atoms
//! that are on different sides.
//!
//! The bacterium carries a genome: a list of genes, each gene being a set of enzymes,
//! each enzyme being one reaction. The chemistry that can happen is
//! environment reactions (outside only) + every enzyme of every gene.
//!
//! Examples of enzymes:
//! * `A1 -> a1`               importer (takes an A1 in, with its molecule)
//! * `a1 -> A1`               exporter
//! * `A1 + b2 -> a2.b3`       receptor b2 takes A1 in, binds it and becomes active b3
//! * `a2.b3 -> a2 + b2`       receptor releases its cargo and goes back to b2

use crate::atoms::is_outside_form;
use crate::molecules::Molecule;
use crate::reactions::{parse_atom_token, Reaction};
use crate::reactor::Reactor;
use rand::prelude::StdRng;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compartment {
    Outside,
    Inside,
}

/// Compartment of an atom key (form << 8 | state).
pub fn compartment(key: u16) -> Compartment {
    if is_outside_form((key >> 8) as u8) { Compartment::Outside } else { Compartment::Inside }
}

/// What an enzyme does with respect to the membrane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnzymeKind {
    /// every atom stays inside
    Internal,
    /// every atom stays outside (not allowed in a genome for now)
    External,
    /// at least one atom goes in, none goes out
    Import,
    /// at least one atom goes out, none goes in
    Export,
    /// one atom goes in while another goes out
    Antiport,
    /// atoms on both sides but nobody crosses (e.g. A1 + b2 -> A2 + b3, a signal)
    Membrane,
}

/// Classify a reaction by the compartments of its atoms.
pub fn classify(r: &Reaction) -> EnzymeKind {
    // (substrate, product) pairs, slot by slot
    let pairs: Vec<(Compartment, Compartment)> = if r.is_production() || r.is_degradation() {
        let c = compartment(r.get_substrate_a()); // both store their atom in slot a
        vec![(c, c)]
    } else if r.is_transform() {
        vec![(compartment(r.get_substrate_a()), compartment(r.get_product_a()))]
    } else {
        vec![
            (compartment(r.get_substrate_a()), compartment(r.get_product_a())),
            (compartment(r.get_substrate_b()), compartment(r.get_product_b())),
        ]
    };

    use Compartment::*;
    let goes_in = pairs.iter().any(|&p| p == (Outside, Inside));
    let goes_out = pairs.iter().any(|&p| p == (Inside, Outside));
    match (goes_in, goes_out) {
        (true, true) => EnzymeKind::Antiport,
        (true, false) => EnzymeKind::Import,
        (false, true) => EnzymeKind::Export,
        (false, false) => {
            if pairs.iter().all(|&(s, _)| s == Inside) {
                EnzymeKind::Internal
            } else if pairs.iter().all(|&(s, _)| s == Outside) {
                EnzymeKind::External
            } else {
                EnzymeKind::Membrane
            }
        }
    }
}

/// One enzyme = one reaction.
#[derive(Debug, Clone)]
pub struct Enzyme {
    pub reaction: Reaction,
    pub kind: EnzymeKind,
}

impl Enzyme {
    pub fn new(reaction: Reaction) -> Result<Enzyme, String> {
        let kind = classify(&reaction);
        if kind == EnzymeKind::External {
            return Err(format!(
                "enzyme '{}' only touches the outside; put it in the environment instead",
                reaction.to_string()
            ));
        }
        Ok(Enzyme { reaction, kind })
    }

    /// Build from the text notation, e.g. `Enzyme::from_rule("A1 -> a1", 0.1)`.
    pub fn from_rule(rule: &str, speed: f64) -> Result<Enzyme, String> {
        Enzyme::new(Reaction::from_rule(0, rule, speed)?)
    }
}

/// A gene is a named set of enzymes.
#[derive(Debug, Clone)]
pub struct Gene {
    pub name: String,
    pub enzymes: Vec<Enzyme>,
}

impl Gene {
    pub fn new(name: &str) -> Gene {
        Gene { name: name.to_string(), enzymes: Vec::new() }
    }

    /// Builder: `Gene::new("receptor").with("A1 + b2 -> a2.b3", 0.01)?.with(...)?`
    pub fn with(mut self, rule: &str, speed: f64) -> Result<Gene, String> {
        self.enzymes.push(Enzyme::from_rule(rule, speed)?);
        Ok(self)
    }
}

#[derive(Debug, Clone, Default)]
pub struct Bacterium {
    pub genes: Vec<Gene>,
}

impl Bacterium {
    pub fn new() -> Bacterium {
        Bacterium { genes: Vec::new() }
    }
    pub fn add_gene(&mut self, gene: Gene) {
        self.genes.push(gene);
    }
}

/// Where a reaction of the reactor comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Origin {
    Environment,
    Gene { gene: usize, enzyme: usize },
}

/// A reactor holding one bacterium in a medium.
///
/// It wraps the standard Gillespie `Reactor`: all atoms (inside and outside) live in the
/// same reactor, the compartment is read from the case of the letter. The reactions are the
/// environment ones plus the enzymes of the bacterium.
pub struct BacteriumReactor {
    reactor: Reactor,
    bacterium: Bacterium,
    /// parallel to `reactor.reactions()`
    origins: Vec<Origin>,
}

impl BacteriumReactor {
    /// Create the reactor and install every enzyme of the bacterium.
    pub fn new(bacterium: Bacterium) -> BacteriumReactor {
        let mut br = BacteriumReactor { reactor: Reactor::new(), bacterium, origins: Vec::new() };
        for (g, gene) in br.bacterium.genes.iter().enumerate() {
            for (e, enzyme) in gene.enzymes.iter().enumerate() {
                br.reactor.add_reaction(enzyme.reaction.clone());
                br.origins.push(Origin::Gene { gene: g, enzyme: e });
            }
        }
        br
    }

    /// Add a reaction of the medium (food supply, dilution, outside chemistry).
    /// It must only involve outside (uppercase) atoms.
    pub fn add_environment_rule(&mut self, rule: &str, speed: f64) -> Result<(), String> {
        let r = Reaction::from_rule(0, rule, speed)?;
        if classify(&r) != EnzymeKind::External {
            return Err(format!("environment reaction '{rule}' must only involve outside (uppercase) atoms"));
        }
        self.reactor.add_reaction(r);
        self.origins.push(Origin::Environment);
        Ok(())
    }

    /// Add `n` free atoms, e.g. `add_atoms("A1", 100)` (outside) or `add_atoms("b2", 10)` (inside).
    pub fn add_atoms(&mut self, atom: &str, n: usize) -> Result<(), String> {
        let key = parse_atom_token(atom.trim())?;
        self.reactor.add_atom_from_keys(key, n);
        Ok(())
    }

    /// Add `n` copies of a linear molecule, e.g. `add_chains("B1.A1", 10)` or `add_chains("c1.c2.c3", 5)`.
    /// All atoms of a molecule must be on the same side.
    pub fn add_chains(&mut self, chain: &str, n: usize) -> Result<(), String> {
        let keys: Vec<u16> = chain.split('.').map(|t| parse_atom_token(t.trim())).collect::<Result<_, _>>()?;
        if keys.iter().any(|&k| compartment(k) != compartment(keys[0])) {
            return Err(format!("molecule '{chain}' must be on one side of the membrane"));
        }
        let mut mol = Molecule::new();
        for (i, &k) in keys.iter().enumerate() {
            mol.add_atom(i as u32, ((k >> 8) as u8, (k & 0xff) as u8));
        }
        for i in 1..keys.len() {
            mol.links.insert(i as u32, (i as u32 - 1, i as u32));
        }
        for _ in 0..n {
            self.reactor.add_molecule(&mol);
        }
        Ok(())
    }

    /// Run Gillespie until `tmax` (can be called repeatedly with increasing `tmax`).
    pub fn run(&mut self, tmax: f64, rng: &mut StdRng) -> f64 {
        self.reactor.main_loop(tmax, rng)
    }

    pub fn time(&self) -> f64 { self.reactor.current_time() }
    pub fn reactor(&self) -> &Reactor { &self.reactor }
    pub fn reactor_mut(&mut self) -> &mut Reactor { &mut self.reactor }
    pub fn bacterium(&self) -> &Bacterium { &self.bacterium }
    pub fn origins(&self) -> &[Origin] { &self.origins }

    /// Number of atoms in a compartment.
    pub fn count(&self, c: Compartment) -> usize {
        self.reactor.iter().filter(|a| compartment(a.get_key()) == c).count()
    }

    /// Atom counts per symbol ("a1", "B3", ...) for one compartment, sorted by symbol.
    pub fn census(&self, c: Compartment) -> BTreeMap<String, usize> {
        let mut m = BTreeMap::new();
        for a in self.reactor.iter().filter(|a| compartment(a.get_key()) == c) {
            *m.entry(a.to_string()).or_insert(0) += 1;
        }
        m
    }

    /// How many times each gene's enzymes fired (indexed like `bacterium().genes`).
    pub fn gene_activity(&self) -> Vec<u64> {
        let mut v = vec![0u64; self.bacterium.genes.len()];
        for (o, &n) in self.origins.iter().zip(self.reactor.fire_counts()) {
            if let Origin::Gene { gene, .. } = o {
                v[*gene] += n;
            }
        }
        v
    }

    /// Short text report: compartments, then genes with their enzymes and fire counts.
    pub fn report(&self) -> String {
        let mut s = format!("t = {:.3}\n", self.time());
        for (label, c) in [("outside", Compartment::Outside), ("inside", Compartment::Inside)] {
            let census: Vec<String> = self.census(c).iter().map(|(k, v)| format!("{k}:{v}")).collect();
            s += &format!("  {label:8} {:5}  [{}]\n", self.count(c), census.join(" "));
        }
        s += &format!("  bonds    {:5}\n", self.reactor.get_link_number());
        let fires = self.reactor.fire_counts();
        for (g, gene) in self.bacterium.genes.iter().enumerate() {
            s += &format!("  gene {} '{}'\n", g, gene.name);
            for (e, enzyme) in gene.enzymes.iter().enumerate() {
                let idx = self.origins.iter()
                    .position(|o| *o == Origin::Gene { gene: g, enzyme: e }).unwrap();
                s += &format!("    {:22} {:9?} fired {}\n",
                    enzyme.reaction.to_string(), enzyme.kind, fires[idx]);
            }
        }
        s
    }
}
