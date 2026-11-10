use rand::{SeedableRng};
use ssa_model::{Reaction};
use ssa_model::reactor::Reactor;
#[cfg(test)]
mod tests {
    use super::*;

    use rand::prelude::StdRng;
    use ssa_model::molecules::Molecule;
    use ssa_model::trace;

    #[test]
    fn test_reactor_create() {
        let seed = 42u64;
        let rng = StdRng::seed_from_u64(seed);

        let mut reactor = Reactor::new();
        let n = 20;
        let form = 0x01u8;
        let state = 0x01u8;
        let key = (form as u16) << 8 | state as u16;
        reactor.add_atom_from_keys(key,n);

        // filter by exact value
        let x = reactor.get_atom_number_from_key(key);
        assert_eq!(x, n);
        // filter with form mask (here 0xff)
        let form = 0xffu8;
        let key = (form as u16) << 8 | state as u16;
        let x = reactor.get_atom_number_from_key(key);
        assert_eq!(x, 0);
        // filter with state mask (should yield 0)
        let state = 0xffu8;
        let form = 0x01u8;
        let key = (form as u16) << 8 | state as u16;
        let x = reactor.get_atom_number_from_key(key);
        assert_eq!(x, 0);


    }

    #[test]
    fn test_reactor_reaction1() {
        let seed = 42u64;
        let mut rng = StdRng::seed_from_u64(seed);

        let mut reactor = Reactor::new();
        let n = 20;
        let form = 0x01u8;
        let state = 0x01u8;
        let key1 = (form as u16) << 8 | state as u16;
        reactor.add_atom_from_keys(key1,n);

        let n = 20;
        let form = 0x01u8;
        let state = 0x02u8;
        let key2 = (form as u16) << 8 | state as u16;

        reactor.add_atom_from_keys(key2,n);
        let reactions = Reaction::news_from_json("data/test_r1.json",1);

        let product_key = reactions[0].product_key().clone();
        let key3 = ((product_key & 0xffff0000) >> 16) as u16;
        let key4 = ((product_key & 0x0000ffff) >> 0) as u16;
        // first initialization tests
        let rs = reactions[0].reaction_speed().clone();
        assert_eq!(rs, 1.0);

        for r in reactions {
            reactor.add_reaction(r);
        }
        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 1);
        assert_eq!(propensity[0], ((n * n) as f64) * rs);
        let a0 = reactor.compute_a0(&propensity);
        assert_eq!(a0, ((n * n) as f64) * rs);

        // apply the only reaction
        reactor.apply_reaction(0, &mut rng);

        reactor.display_stats();

        let x1 = reactor.get_atom_number_from_key(key1);
        assert_eq!(x1, n-1);

        let x2 = reactor.get_atom_number_from_key(key2);
        assert_eq!(x2, n-1);

        let x3 = reactor.get_atom_number_from_key(key3);
        assert_eq!(x3, 1);

        let x4 = reactor.get_atom_number_from_key(key4);
        assert_eq!(x4, 1);

        let x = reactor.iter_links().count();
        assert_eq!(x, 1);
    }
    #[test]
    fn test_reactor_reaction2() {
        let seed = 42u64;
        let mut rng = StdRng::seed_from_u64(seed);
        let mut reactor = Reactor::new();
        let mut molecule = Molecule::new();
        molecule.load_from_file("data/reactor1.mol");
        reactor.add_molecule(&molecule);
        let reactions = Reaction::news_from_json("data/test_r2.json",1);
        for r in reactions {
            reactor.add_reaction(r);
        }

        let key = 0x01010102;
        assert_eq!(reactor.get_link_from_key(key).len(), 1);

        let key = 0x01010103;
        assert_eq!(reactor.get_link_from_key(key).len(), 1);

        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 1);
        assert_eq!(propensity[0], 1.0f64);
        reactor.apply_reaction(0, &mut rng);
        reactor.display_stats();
        let key1 = 0x0109;
        let key2 = 0x010A;
        let key3 = 0x0101;
        let key4 = 0x0102;
        assert_eq!(reactor.get_atom_number_from_key(key1), 1);
        assert_eq!(reactor.get_atom_number_from_key(key2), 1);
        assert_eq!(reactor.get_atom_number_from_key(key3), 1);
        assert_eq!(reactor.get_atom_number_from_key(key4), 0);

        for x in reactor.iter_links() {
            trace!(x.key());
        }

        let key = 0x01010102;
        assert_eq!(reactor.get_link_from_key(key).len(), 0);
        let key = 0x01010103;
        assert_eq!(reactor.get_link_from_key(key).len(), 1);
    }

    #[test]
    fn test_reactor_reaction3() {
        let seed = 42u64;
        let mut rng = StdRng::seed_from_u64(seed);
        let mut reactor = Reactor::new();
        let mut molecule = Molecule::new();
        molecule.load_from_file("data/reactor2.mol");
        reactor.add_molecule(&molecule);
        let reactions = Reaction::news_from_json("data/test_r3.json",1);
        for r in reactions {
            reactor.add_reaction(r);
        }



        let key = 0xFF01;
        assert_eq!(reactor.get_atom_number_from_key(key), 0);
        let key = 0xFF02;
        assert_eq!(reactor.get_atom_number_from_key(key), 0);

        let key = 0x0402;
        assert_eq!(reactor.get_atom_number_from_key(key), 1);

        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 256);
        //assert_eq!(propensity[0], 4.0f64);
        let a0 = reactor.compute_a0(&propensity);
        let t = reactor.compute_next_reaction_time(&propensity, & mut rng);
        let id = reactor.compute_next_reaction(&propensity, & mut rng).unwrap();
        reactor.apply_reaction(id, &mut rng);
        reactor.display_stats();
        let key1 = 0xFF01;
        let key2 = 0xFF02;
        let key3 = 0xFF09;
        let key4 = 0xFF0A;
        assert_eq!(reactor.get_atom_number_from_key(key1), 2);
        assert_eq!(reactor.get_atom_number_from_key(key2), 3);
        assert_eq!(reactor.get_atom_number_from_key(key3), 1);
        assert_eq!(reactor.get_atom_number_from_key(key4), 1);

        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 256);
        //assert_eq!(propensity[0], 2.0f64);
    }
    #[test]
    fn test_reactor_reaction4() {
        let seed = 42u64;
        let mut rng = StdRng::seed_from_u64(seed);
        let mut reactor = Reactor::new();
        let mut molecule = Molecule::new();
        molecule.load_from_file("data/reactor3.mol");
        reactor.add_molecule(&molecule);
        let reactions = Reaction::news_from_json("data/test_r4.json",1);
        for r in reactions {
            reactor.add_reaction(r);
        }

        let key = 0xEF01;
        assert_eq!(reactor.get_atom_number_from_key(key), 3);
        let key = 0xEF02;
        assert_eq!(reactor.get_atom_number_from_key(key), 4);

        let key = 0x0402;
        assert_eq!(reactor.get_atom_number_from_key(key), 1);

        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 1);
        assert_eq!(propensity[0], 12.0f64);

        reactor.apply_reaction(0, &mut rng);
        reactor.display_stats();
        let key1 = 0xFF01;
        let key2 = 0xEF02;
        let key3 = 0xFF09;
        let key4 = 0xFF0A;
        assert_eq!(reactor.get_atom_number_from_key(key1), 2);
        assert_eq!(reactor.get_atom_number_from_key(key2), 3);
        assert_eq!(reactor.get_atom_number_from_key(key3), 1);
        assert_eq!(reactor.get_atom_number_from_key(key4), 1);

        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 1);
        let a = (reactor.get_atom_number_from_key(key1) * reactor.get_atom_number_from_key(key2)) as f64;
        assert_eq!(propensity[0], a);
    }
    #[test]
    fn test_reactor_propagation1() {
        let seed = 42u64;
        let mut rng = StdRng::seed_from_u64(seed);
        let mut reactor = Reactor::new();
        let mut molecule = Molecule::new();
        molecule.load_from_file("data/reactor4.mol");
        reactor.add_molecule(&molecule);



        let rv = Reaction::news_from_json("data/test_r5.json",1);
        for r in rv {
            reactor.add_reaction(r);
        }
        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 3);
        assert_eq!(propensity[0], 1f64);
        for x in &propensity {println!("{:}",*x);}
        reactor.apply_reaction(0, &mut rng);
        reactor.display_stats();
        for l in reactor.iter_links() {println!("{:} {:}",l.links[0], l.links[1]);}
        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 3);
        assert_eq!(propensity[1], 1f64);
        for x in &propensity {println!("{:}",*x);}
        reactor.apply_reaction(1, &mut rng);
        reactor.display_stats();
        for l in reactor.iter_links() {println!("{:} {:}",l.links[0], l.links[1]);}
        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 3);
        for x in &propensity {println!("{:}",*x);}
        assert_eq!(propensity[2], 1f64);
        reactor.apply_reaction(2, &mut rng);
        reactor.display_stats();
        for l in reactor.iter_links() {println!("{:} {:}",l.links[0], l.links[1]);}
    }

    #[test]
    fn test_reactor_propagation2() {
        let seed = 42u64;
        let mut rng = StdRng::seed_from_u64(seed);
        let mut reactor = Reactor::new();
        let mut molecule = Molecule::new();
        molecule.load_from_file("data/reactor5.mol");
        reactor.add_molecule(&molecule);
        let reactions = Reaction::news_from_json("data/test_r6.json",1);
        for r in reactions {
            println!("{:}",r.to_string());
            reactor.add_reaction(r);
        }
        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 16 * 3);
        let a0 = reactor.compute_a0(&propensity);
        let t = reactor.compute_next_reaction_time(&propensity, & mut rng);
        let id = reactor.compute_next_reaction(&propensity, & mut rng).unwrap();
        reactor.apply_reaction(id, &mut rng);
        reactor.display_stats();

        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 16 * 3);
        let a0 = reactor.compute_a0(&propensity);
        let t = reactor.compute_next_reaction_time(&propensity, & mut rng);
        let id = reactor.compute_next_reaction(&propensity, & mut rng).unwrap();
        reactor.apply_reaction(id, &mut rng);
        reactor.display_stats();

        let propensity = reactor.get_propensity_vector();
        assert_eq!(propensity.len(), 16 * 3);
        let a0 = reactor.compute_a0(&propensity);
        let t = reactor.compute_next_reaction_time(&propensity, & mut rng);
        let id = reactor.compute_next_reaction(&propensity, & mut rng).unwrap();
        reactor.apply_reaction(id, &mut rng);
        reactor.display_stats();
    }
}