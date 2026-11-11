use rand::SeedableRng;
use rand::prelude::StdRng;
use ssa_model::bacterium::{classify, Bacterium, BacteriumReactor, Compartment, Enzyme, EnzymeKind, Gene};
use ssa_model::reactions::Reaction;

fn run_to_end(br: &mut BacteriumReactor, seed: u64) {
    let mut rng = StdRng::seed_from_u64(seed);
    br.run(1e6, &mut rng); // stops by itself once no reaction is possible
}

#[test]
fn parse_rules_roundtrip() {
    for (rule, printed) in [
        ("A1 -> a1", "A1->a1"),
        ("a1 -> A1", "a1->A1"),
        ("A1 + b2 -> a2.b3", "A1+b2->a2.b3"),
        ("a2.b3 -> a2 + b2", "a2.b3->a2+b2"),
        ("0 -> A1", "0->A1"),
        ("A1 -> 0", "A1->0"),
        ("a1 + b2 -> a3.b4", "a1+b2->a3.b4"), // old notation unchanged
    ] {
        let r = Reaction::from_rule(0, rule, 1.0).unwrap();
        assert_eq!(r.to_string(), printed, "rule {rule}");
    }
    assert!(Reaction::from_rule(0, "a1 -> b1", 1.0).is_err(), "form change must be rejected");
    assert!(Reaction::from_rule(0, "A1 -> b1", 1.0).is_err(), "form change must be rejected");
    assert!(Reaction::from_rule(0, "a1 -> a2 + b3", 1.0).is_err(), "arity mismatch must be rejected");
    assert!(Reaction::from_rule(0, "A1 + b2 -> A2.b3", 1.0).is_err(), "bond across membrane");
    assert!(Reaction::from_rule(0, "A1.b2 -> a1 + b2", 1.0).is_err(), "bond across membrane");
}

#[test]
fn classify_enzymes() {
    let k = |rule: &str| classify(&Reaction::from_rule(0, rule, 1.0).unwrap());
    assert_eq!(k("A1 -> a1"), EnzymeKind::Import);
    assert_eq!(k("a1 -> A1"), EnzymeKind::Export);
    assert_eq!(k("A1 + b2 -> a2.b3"), EnzymeKind::Import);
    assert_eq!(k("a2.b3 -> a2 + b2"), EnzymeKind::Internal);
    assert_eq!(k("A1 + b2 -> A2 + b3"), EnzymeKind::Membrane);
    assert_eq!(k("a1 + B2 -> A1 + b2"), EnzymeKind::Antiport);
    assert_eq!(k("A1 + B2 -> A3 + B4"), EnzymeKind::External);
    assert!(Enzyme::from_rule("A1 + B2 -> A3 + B4", 1.0).is_err());
}

#[test]
fn importer_takes_free_atoms_in() {
    let mut b = Bacterium::new();
    b.add_gene(Gene::new("importer").with("A1 -> a1", 1.0).unwrap());
    let mut br = BacteriumReactor::new(b);
    br.add_atoms("A1", 50).unwrap();
    run_to_end(&mut br, 1);
    assert_eq!(br.count(Compartment::Outside), 0);
    assert_eq!(br.census(Compartment::Inside).get("a1"), Some(&50));
    assert_eq!(br.gene_activity(), vec![50]);
}

#[test]
fn whole_molecule_crosses() {
    // B1.A1 outside, only A1 is recognised by the importer: B1 must follow and keep its state
    let mut b = Bacterium::new();
    b.add_gene(Gene::new("importer").with("A1 -> a1", 1.0).unwrap());
    let mut br = BacteriumReactor::new(b);
    br.add_chains("B1.A1", 10).unwrap();
    br.add_chains("C5.A1.D7", 5).unwrap();
    run_to_end(&mut br, 2);
    let inside = br.census(Compartment::Inside);
    assert_eq!(br.count(Compartment::Outside), 0);
    assert_eq!(inside.get("a1"), Some(&15));
    assert_eq!(inside.get("b1"), Some(&10));
    assert_eq!(inside.get("c5"), Some(&5));
    assert_eq!(inside.get("d7"), Some(&5));
    assert_eq!(br.reactor().get_link_number(), 20);
}

#[test]
fn whole_molecule_goes_out() {
    let mut b = Bacterium::new();
    b.add_gene(Gene::new("exporter").with("a3 -> A3", 1.0).unwrap());
    let mut br = BacteriumReactor::new(b);
    br.add_chains("b1.b2.a3.b4", 4).unwrap();
    br.add_chains("c1.c2", 3).unwrap(); // not connected to any a3: stays in
    run_to_end(&mut br, 3);
    let outside = br.census(Compartment::Outside);
    assert_eq!(outside.get("A3"), Some(&4));
    assert_eq!(outside.get("B1"), Some(&4));
    assert_eq!(outside.get("B4"), Some(&4));
    assert_eq!(br.count(Compartment::Inside), 6);
}

#[test]
fn receptor_imports_and_releases() {
    let mut b = Bacterium::new();
    b.add_gene(
        Gene::new("receptor")
            .with("A1 + b2 -> a2.b3", 0.01).unwrap()
            .with("a2.b3 -> a2 + b2", 1.0).unwrap(),
    );
    let mut br = BacteriumReactor::new(b);
    br.add_atoms("A1", 100).unwrap();
    br.add_atoms("b2", 10).unwrap();
    run_to_end(&mut br, 7);
    let inside = br.census(Compartment::Inside);
    assert_eq!(br.count(Compartment::Outside), 0);
    assert_eq!(inside.get("a2"), Some(&100));
    assert_eq!(inside.get("b2"), Some(&10));
    assert_eq!(br.reactor().get_link_number(), 0);
}

#[test]
fn receptor_drags_cargo_molecule() {
    // the A1 grabbed by the receptor is bonded to C9 outside: C9 comes in as c9
    let mut b = Bacterium::new();
    b.add_gene(
        Gene::new("receptor")
            .with("A1 + b2 -> a2.b3", 1.0).unwrap()
            .with("a2.b3 -> a2 + b2", 1.0).unwrap(),
    );
    let mut br = BacteriumReactor::new(b);
    br.add_chains("C9.A1", 20).unwrap();
    br.add_atoms("b2", 3).unwrap();
    run_to_end(&mut br, 11);
    let inside = br.census(Compartment::Inside);
    assert_eq!(br.count(Compartment::Outside), 0);
    assert_eq!(inside.get("c9"), Some(&20));
    assert_eq!(inside.get("a2"), Some(&20));
    assert_eq!(inside.get("b2"), Some(&3));
    assert_eq!(br.reactor().get_link_number(), 20); // only the C9.A1 bonds remain
}

#[test]
fn bond_rule_written_in_reverse_order() {
    // bond created as (a, b) but the release rule is written b.a
    let mut b = Bacterium::new();
    b.add_gene(
        Gene::new("receptor")
            .with("A1 + b2 -> a2.b3", 1.0).unwrap()
            .with("b3.a2 -> b2 + a5", 1.0).unwrap(),
    );
    let mut br = BacteriumReactor::new(b);
    br.add_atoms("A1", 20).unwrap();
    br.add_atoms("b2", 5).unwrap();
    run_to_end(&mut br, 3);
    let inside = br.census(Compartment::Inside);
    assert_eq!(inside.get("a5"), Some(&20));
    assert_eq!(inside.get("b2"), Some(&5));
}

#[test]
fn environment_rules_must_be_outside() {
    let mut br = BacteriumReactor::new(Bacterium::new());
    assert!(br.add_environment_rule("0 -> A1", 1.0).is_ok());
    assert!(br.add_environment_rule("A1 -> 0", 0.1).is_ok());
    assert!(br.add_environment_rule("0 -> a1", 1.0).is_err());
    assert!(br.add_environment_rule("A1 -> a1", 1.0).is_err());
    assert!(br.add_chains("B1.a1", 1).is_err());
}
