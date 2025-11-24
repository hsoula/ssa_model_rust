pub mod atoms;
pub mod reactions;
mod reactor;
mod links;
mod utils;

pub use atoms::Atom;
pub use reactions::Reaction;

#[cfg(test)]
mod tests {
    use crate::{Atom, Reaction};

    #[test]
    fn test_atom_create_id() {
        let id : u32 = 1;
        let form : u8 = 0x0f;
        let state : u8 = 0x01;
        let atom = Atom::new(id, form, state);
        assert_eq!(atom.id(), id);
    }
    #[test]
    fn test_atom_create_form() {
        let id : u32 = 1;
        let form : u8 = 0x0f;
        let state : u8 = 0x01;
        let atom = Atom::new(id, form, state);
        assert_eq!(atom.form(), form);
    }
    #[test]
    fn test_atom_create_state() {
        let id : u32 = 1;
        let form : u8 = 0x0f;
        let state : u8 = 0x01;
        let atom = Atom::new(id, form, state);
        assert_eq!(atom.state(), state);
    }
    #[test]
    fn test_atom_create_key() {
        let id : u32 = 1;
        let form : u8 = 0x0f;
        let state : u8 = 0x01;
        let atom = Atom::new(id, form, state);
        assert_eq!(atom.get_key(), 0x0f01u16);
    }
    #[test]
    fn test_atom_insert(){
        let id : u32 = 1;
        let form : u8 = 0x0f;
        let state : u8 = 0x01;
        let mut atom = Atom::new(id, form, state);
        atom.link(10);
        assert_eq!(atom.get_nb_link(), 1);
        assert_eq!(atom.is_link(10), true);
        assert_eq!(atom.is_link(1), false);
        atom.link(10);
        assert_eq!(atom.get_nb_link(), 1);
        atom.link(20);
        assert_eq!(atom.get_nb_link(), 2);

    }
    #[test]
    fn test_atom_unlink(){
        let id : u32 = 1;
        let form : u8 = 0x0f;
        let state : u8 = 0x01;
        let mut atom = Atom::new(id, form, state);
        atom.link(10);
        assert_eq!(atom.get_nb_link(), 1);
        assert_eq!(atom.is_link(10), true);
        atom.link(20);
        assert_eq!(atom.is_link(20), true);
        assert_eq!(atom.get_nb_link(), 2);
        atom.unlink(20);
        assert_eq!(atom.is_link(20), false);
        assert_eq!(atom.get_nb_link(), 1);
    }


    #[test]
    fn test_reaction_creation(){
        let id : u32 = 1;
        let r = Reaction::new(id);
        assert_eq!(r.is_product_contact(), false);
    }
    #[test]
    fn test_reaction_creation_json() {
        let r = Reaction::new_from_json("data/r1.json",1);
        //assert_eq!(r.reaction_type(),):
        println!("{:?}", r.reaction_type());
        assert_eq!(r.is_product_contact(), false);
        assert_eq!(r.is_substrate_contact(), true);
        assert_eq!(r.get_substrate_a(), 0x0100);
        assert_eq!(r.get_substrate_b(), 0x0101);
    }
}
