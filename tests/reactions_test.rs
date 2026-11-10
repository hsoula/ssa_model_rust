
use ssa_model::{Atom, Reaction, };
use ssa_model::reactor::Reactor;
use ssa_model::reactions::ReactionFlag;
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::NamedTempFile;
    use ssa_model::trace;

    #[test]
    fn test_reaction_creation_json() {
        let reactions = Reaction::news_from_json("data/r1.json", 1);
        for r in reactions {
            //assert_eq!(r.reaction_type(),):
            println!("{:?}", r.reaction_type());
            assert_eq!(r.is_product_contact(), false);
            assert_eq!(r.is_substrate_contact(), true);
            assert_eq!(r.get_substrate_a(), 0x0101);
            assert_eq!(r.get_substrate_b(), 0x0102);
        }
    }

    #[test]
    fn test_reactor_creation() {
        let mut reactor = Reactor::new();
        let n = 100;
        let form = 0x02u8;
        let state = 0x01u8;
        let key = (form as u16) << 8 | state as u16;
        assert_eq!(key, 0x0201u16);
        let key2 = (state as u16) << 8 | form as u16;
        let a = Atom::new(form, state);
        assert_eq!(a.get_key(), key);
        reactor.add_atom_from_keys(key, n);
        assert_eq!(reactor.get_atom_number_from_key(key), n);
        assert_eq!(reactor.get_atom_number_from_key(key2), 0);
        reactor.add_atom_from_keys(key2, n);
        assert_eq!(reactor.get_atom_number_from_key(key2),n);
        assert_eq!(reactor.get_atom_number() ,2 * n);
    }

    // Test data helper

    #[test]
    fn test_reaction_new() {
        let reaction = Reaction::new(42);

        assert_eq!(reaction.id(), 42);
        assert_eq!(reaction.reaction_type(), 0);
        assert_eq!(reaction.reaction_speed(), 1.0);
        assert_eq!(reaction.substrate_key(), 0);
        assert_eq!(reaction.product_key(), 0);
    }

    #[test]
    fn test_reaction_new_from_json() {
        let reactions = Reaction::news_from_json("data/test_r0.json", 1);
        for reaction in reactions {
            assert_eq!(reaction.id(), 1);
            assert_eq!(reaction.reaction_type(), 3);
            assert_eq!(reaction.reaction_speed(), 2.5);
            assert_eq!(reaction.substrate_key(), 0x41420000 | 0x4344);
            assert_eq!(reaction.product_key(), 0x45460000 | 0x4748);
        }
    }


    #[test]
    fn test_is_substrate_contact() {
        let mut r = Reaction::new(1);

        // Test when flag is set
        r = r.set_reaction_type(ReactionFlag::SubstrateContact as u8);
        assert!(r.is_substrate_contact());

        // Test when flag is not set
        r = r.set_reaction_type( 0);
        assert!(!r.is_substrate_contact());

        // Test when other flags are set but not substrate contact
        r = r.set_reaction_type(ReactionFlag::ProductContact as u8);
        assert!(!r.is_substrate_contact());

        // Test when both flags are set
        r = r.set_reaction_type((ReactionFlag::SubstrateContact as u8) | (ReactionFlag::ProductContact as u8));
        assert!(r.is_substrate_contact());
    }

    #[test]
    fn test_is_product_contact() {
        let mut r = Reaction::new(1);

        // Test when flag is set
        r = r.set_reaction_type(ReactionFlag::ProductContact as u8);
        assert!(r.is_product_contact());

        // Test when flag is not set
        r =  r.set_reaction_type(0);
        assert!(!r.is_product_contact());

        // Test when other flags are set but not product contact
        r=r.set_reaction_type( ReactionFlag::SubstrateContact as u8);
        assert!(!r.is_product_contact());

        // Test when both flags are set
        r=r.set_reaction_type((ReactionFlag::SubstrateContact as u8) | (ReactionFlag::ProductContact as u8));
        assert!(r.is_product_contact());
    }


    #[test]
    fn test_get_substrate_a() {
        let mut r = Reaction::new(1);
        r = r.set_substrate_key(0x12345678);

        assert_eq!(r.get_substrate_a(), 0x1234);
    }

    #[test]
    fn test_get_substrate_b() {
        let mut r = Reaction::new(1);
        r = r.set_substrate_key(0x12345678);

        assert_eq!(r.get_substrate_b(), 0x5678);
    }

    #[test]
    fn test_get_product_a() {
        let mut r = Reaction::new(1);
        r = r.set_product_key( 0xABCDEF12);

        assert_eq!(r.get_product_a(), 0xABCD);
    }

    #[test]
    fn test_get_product_b() {
        let mut r = Reaction::new(1);
        r = r.set_product_key(0xABCDEF12);

        assert_eq!(r.get_product_b(), 0xEF12);
    }

    #[test]
    fn test_to_string_basic() {
        let mut reaction = Reaction::new(1);
        reaction = reaction.set_substrate_key ((0x0161 << 16) | 0x0262); // s1='a' (97), s2='b' (98)
        reaction = reaction.set_product_key ((0x0363 << 16) | 0x0464);   // p1='c' (99), p2='d' (100)
        reaction = reaction.set_reaction_type (0);

        let result = reaction.to_string();
        assert_eq!(result, "a97.b98->c99.d100");
    }

    #[test]
    fn test_to_string_with_substrate_contact() {
        let mut r = Reaction::new(1);
        r = r.set_substrate_key((0x0161 << 16) | 0x0262);
        r = r.set_product_key((0x0263 << 16) | 0x0364);
        r = r.set_reaction_type(ReactionFlag::SubstrateContact as u8);

        let result = r.to_string();
        assert_eq!(result, "a97+b98->b99.c100");
    }

    #[test]
    fn test_to_string_with_product_contact() {
        let mut r = Reaction::new(1);
        r= r.set_substrate_key((0x0161 << 16) | 0x0262);
        r = r.set_product_key((0x0363 << 16) | 0x0464);
        r= r.set_reaction_type(ReactionFlag::ProductContact as u8);

        let result = r.to_string();
        assert_eq!(result, "a97.b98->c99+d100");
    }

    #[test]
    fn test_to_string_with_both_contacts() {
        let mut r = Reaction::new(1);
        r = r.set_substrate_key((0x0161 << 16) | 0x0262);
        r = r.set_product_key((0x0363 << 16) | 0x0464);
        r = r.set_reaction_type((ReactionFlag::SubstrateContact as u8) | (ReactionFlag::ProductContact as u8));

        let result = r.to_string();
        assert_eq!(result, "a97+b98->c99+d100");
    }

    #[test]
    fn test_to_string_with_high_bytes() {
        let mut r = Reaction::new(1);
        // High bytes > 200 should produce '*'
        r = r.set_substrate_key((0xC961 << 16) | 0xC862); // 0xC8 = 200
        r = r.set_product_key((0xC963 << 16) | 0xCA64);   // 0xC9 = 201, 0xCA = 202
        r = r.set_reaction_type(0);

        let result = r.to_string();
        // Should show '*' for high bytes > 200
        assert_eq!(result, ")97.(98->)99.*100");
    }

    #[test]
    fn test_to_string_edge_cases() {
        let mut r = Reaction::new(1);

        // Test with minimum values
        r = r.set_substrate_key((0xFF01 << 16) | 0xFE01);
        r = r.set_product_key((0xFF01 << 16) | 0xFE01);
        r = r.set_reaction_type(0);

        let result = r.to_string();
        assert_eq!(result, "*1.?1->*1.?1"); // 96 + 0 = '`'

        // Test with maximum normal values (below 200)
        r = r.set_substrate_key((0x01C7 << 16) | 0x01C7); // 199
        r = r.set_product_key((0x01C7 << 16) | 0x01C7);
        let result = r.to_string();
        assert_eq!(result, "a199.a199->a199.a199"); // 96 + 199 = 295, which wraps around
    }

    #[test]
    fn test_reaction_flag_values() {
        // Verify the flag values are as expected
        assert_eq!(ReactionFlag::SubstrateContact as u8, 0b10);
        assert_eq!(ReactionFlag::ProductContact as u8, 0b1);
    }

    #[test]
    fn test_multiple_reactions() {
        let reactions: Vec<Reaction> = (0..3)
            .map(|i| {
                let mut r = Reaction::new(i);
                r = r.set_substrate_key ((i as u32) << 16);
                r = r.set_product_key ((i as u32) << 8);
                r
            })
            .collect();

        for (i, reaction) in reactions.iter().enumerate() {
            assert_eq!(reaction.id(), i as u32);
            trace!(i);
            assert_eq!(reaction.get_substrate_a(), (i as u16) << 0);
            assert_eq!(reaction.get_product_b(), (i as u16) << 8);
        }
    }

    #[test]
    #[should_panic]
    fn test_new_from_json_invalid_file() {
        Reaction::news_from_json("nonexistent_file.json", 1);
    }

    // Test JSON structure validation
    #[test]
    fn test_json_reaction_structure() {

        let reactions = Reaction::news_from_json("data/test_r0.json", 1);
        for reaction in reactions {
            // Verify all fields are correctly parsed
            assert_eq!(reaction.reaction_type(), 3);
            assert_eq!(reaction.reaction_speed(), 2.5);
            assert_eq!(reaction.get_substrate_a(), 0x4142); // 0x4142 high byte
            assert_eq!(reaction.get_substrate_b(), 0x4344);
            assert_eq!(reaction.get_product_a(), 0x4546); // 0x4546 high byte
            assert_eq!(reaction.get_product_b(), 0x4748);
        }
    }
}
