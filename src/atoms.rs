use std::collections::HashSet;
use serde_json::Value::String;
use utils::{ObjectId, generate_id};
use crate::{trace, utils};
///
/// Struct Atom
///
/// # Summary
/// main class of atoms; for ssa algorithm only the number is required but here we need to track
/// individual particle to create "hubryd" atom structure using links
///
/// # Members
///
/// id: unique ID for reactor database
/// form: the type of the atom 255 max
/// state: the state of atoms 255 max
/// link: the hashset of ids of link atoms
///

/// Form bit marking an atom as *outside* the bacterium.
/// Printed as an uppercase letter: form 1 = `a` (inside), form 1|OUTSIDE = `A` (outside).
/// Crossing the membrane flips this bit and leaves the state untouched.
pub const OUTSIDE: u8 = 0x40;

/// Wildcard forms used by the JSON loader (0xEF, 0xFE, 0xFF) are not compartments.
pub fn is_wildcard_form(form: u8) -> bool {
    form >= 0xE0
}

pub fn is_outside_form(form: u8) -> bool {
    !is_wildcard_form(form) && form & OUTSIDE != 0
}

/// Letter for a form: `a`..`z` inside, `A`..`Z` outside, `*` / `?` for wildcards.
pub fn form_symbol(form: u8) -> char {
    match form {
        0xff => '*',
        0xfe => '?',
        f if is_outside_form(f) && (f & !OUTSIDE) >= 1 && (f & !OUTSIDE) <= 26 => ((f & !OUTSIDE) + 64) as char,
        f => f.wrapping_add(96) as char,
    }
}

/// Symbol of an atom key (form << 8 | state), e.g. `a3`, `B12`.
pub fn key_symbol(key: u16) -> std::string::String {
    format!("{}{}", form_symbol((key >> 8) as u8), key & 0xff)
}

pub struct Atom {
    pub id: ObjectId,
    form: u8,
    state: u8,
    link: HashSet<ObjectId>,
}


impl Atom {
    /// Create Atom
    /// # Parameters
    /// * id: unique id fof the reactor
    /// * form: the type of the atom [0,255]
    /// * state: the state of the atom [0, 255]
    pub fn new(form: u8, state: u8) -> Atom {
        let id = generate_id();
        Atom{id:id, form:form, state:state, link:HashSet::new()}
    }

    pub fn new_from_key(key:u16) -> Atom {
        let id = generate_id();
        let form = ((key & 0xff00) >> 8) as u8;
        let state = (key & 0xff) as u8;
        Atom{id:id, form:form, state:state, link:HashSet::new()}
    }

    /// getter for the key
    /// high bit -> form
    /// low bit -> state

    pub fn get_key(&self) -> u16 {
        (self.form as u16) << 8 | self.state as u16
        //u16::from_be_bytes([self.form, self.state])
    }
    /// getter for the id
    pub fn id(&self) -> ObjectId {
        self.id
    }
    /// getter for the form
    pub fn form(&self) -> u8 {
        self.form
    }
    /// getter for the state
    pub fn state(&self) -> u8 {
        self.state
    }
    /// method to link
    /// # argument:
    /// id of atom in reactor database
    /// # warning:
    /// insert only once and prevent self inserting
    ///
    pub fn link(&mut self, object_id: ObjectId)  {
        if object_id != self.id {
            assert!(!self.is_link(object_id));
            self.link.insert(object_id);
        }
    }
    /// method to unlink
    /// # argument:
    /// id of atom in reactor database
    /// # warning:
    /// quiet removing does nothing
    /// # TODO:
    ///  - return the bool is the removed id was in the set for test
    /// mutable action
    pub fn unlink(&mut self, id: ObjectId) {
        self.link.remove(&id);
    }
    /// method to check link
    /// # argument:
    /// id of atom in reactor database
    /// mutable action
    pub fn is_link(&self, id: ObjectId) -> bool {
        self.link.contains(&id)
    }
    /// compute the number of atoms linked
    ///
    pub fn get_nb_link(&self) -> u32 {
        self.link.len() as u32
    }
    /// set state of the atom
    /// # Arguments
    /// * state -> new state
    /// mutable action
    pub fn set_state(&mut self, state: u8) {
        self.state = state;
    }
    /// set form (only used to flip the OUTSIDE bit when crossing the membrane)
    pub fn set_form(&mut self, form: u8) {
        self.form = form;
    }
    pub fn is_outside(&self) -> bool {
        is_outside_form(self.form)
    }

    pub fn get_links(&self) -> &HashSet<ObjectId> {
        &self.link
    }
    pub fn get_linked(&self) -> Vec<ObjectId> {
        self.link.iter().map(|id| (*id).clone()).collect()
    }

    pub fn unlink_all(&mut self)-> Vec<ObjectId> {
        let v = self.link.iter().map(|id| (*id).clone()).collect();
        self.link.clear();
        v
    }
    pub fn to_string(&self) -> std::string::String {
        key_symbol(self.get_key())
    }
}

