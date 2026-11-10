//use crate::trace;
use crate::utils::{ObjectId, generate_id};
use crate::atoms::Atom;

pub struct Link {
    pub id : ObjectId,
    pub links : [ObjectId; 2],
    pub key : u32
}

pub fn get_key_from_atom(a: &Atom, b:&Atom) -> u32 {
    let a = a.get_key();
    let b = b.get_key();
    get_key(a, b)
}

pub fn get_key(a: u16, b:u16) -> u32 {
    let key = ((a as u32) << 16) | b as u32;
    key
}

impl Link {
    pub fn new(a_id : ObjectId, a_key: u16, b_id: ObjectId, b_key: u16) -> Link {
        let id = generate_id();
        let key = get_key(a_key, b_key);
        Link{id: id, links: [a_id, b_id], key: key}
    }
    pub fn update_key(&mut self, a: u16, b:u16) {
        self.key = get_key(a, b);
    }
    pub fn get_id(&self) -> ObjectId {
        self.id
    }
    pub fn id(&self) -> ObjectId {
        self.id
    }
    pub fn key(&self) -> u32 { self.key }
    pub fn get_substrate_key(&self) -> u16 {
        (self.key >> 16) as u16
    }
    pub fn key_rev(&self) -> u32 {
        (self.key & 0xffff0000) >> 16 | (self.key & 0x0000ffff) << 16
    }
}