use std::collections::HashSet;
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

pub struct Atom {
    id: u32,
    form: u8,
    state: u8,
    link: HashSet<u32>,
}


impl Atom {
    /// Create Atom
    /// # Parameters
    /// * id: unique id fof the reactor
    /// * form: the type of the atom [0,255]
    /// * state: the state of the atom [0, 255]
    pub fn new(id: u32, form: u8, state: u8) -> Atom {
        Atom{id:id, form:form, state:state, link:HashSet::new()}
    }
    /// getter for the key
    /// high bit -> form
    /// low bit -> state

    pub fn get_key(&self) -> u16 {
        u16::from_be_bytes([self.form, self.state])
    }
    /// getter for the id
    pub fn id(&self) -> u32 {
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
    pub fn link(&mut self, id: u32)  {
        if id != self.id {
            self.link.insert(id);
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
    pub fn unlink(&mut self, id: u32) {
        self.link.remove(&id);
    }
    /// method to check link
    /// # argument:
    /// id of atom in reactor database
    /// mutable action
    pub fn is_link(&self, id: u32) -> bool {
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
}

