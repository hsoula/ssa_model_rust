

use ssa_model::{Atom};
#[test]
fn test_atom_create_form() {
    let form: u8 = 0x0f;
    let state: u8 = 0x01;
    let atom = Atom::new(form, state);
    assert_eq!(atom.form(), form);
}
#[test]
fn test_atom_create_state() {
    let form: u8 = 0x0f;
    let state: u8 = 0x01;
    let atom = Atom::new(form, state);
    assert_eq!(atom.state(), state);
}
#[test]
fn test_atom_create_key() {
    let form: u8 = 0x0f;
    let state: u8 = 0x01;
    let atom = Atom::new(form, state);
    assert_eq!(atom.get_key(), 0x0f01u16);
}
#[test]
fn test_atom_create_key2() {
    let form: u8 = 0x0f;
    let state: u8 = 0x01;
    let key = (form as u16) << 8 | state as u16;
    println!("{:4x}", key);
    let atom = Atom::new_from_key(key);
    println!("{:4x}", atom.get_key());
    assert_eq!(atom.get_key(), 0x0f01u16);
    assert_eq!(atom.form(), form);
    assert_eq!(atom.state(), state);
}
#[test]
fn test_atom_insert() {
    let form: u8 = 0x0f;
    let state: u8 = 0x01;
    let mut atom = Atom::new(form, state);
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
fn test_atom_unlink() {
    let form: u8 = 0x0f;
    let state: u8 = 0x01;
    let mut atom = Atom::new(form, state);
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
