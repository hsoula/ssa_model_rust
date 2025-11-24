use crate::atoms::Atom;
use crate::links::Link;
use crate::reactions::Reaction;
use std::collections::HashMap;
use rand::Rng;
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;

pub struct Reactor {
    nb_current_atoms: usize,
    links: Vec<Link>,
    atoms: Vec<Atom>,
    reactions: Vec<Reaction>,
    current_time : f64,
    key_count_cache: HashMap<u16, usize>,
    link_count_cache: HashMap<u32, usize>,
    link_remove : Vec<u32>
}

impl Reactor {
    pub fn new() -> Reactor {
        Reactor{nb_current_atoms: 0,
            links: Vec::new(),
            atoms: Vec::new(),
            reactions: Vec::new(),
            current_time: 0f64,
            key_count_cache: HashMap::new(),
            link_count_cache: HashMap::new(),
            link_remove: Vec::new()}
    }
    pub fn get_collision_reactions(&self) -> Vec<&Reaction> {
        self.reactions.iter()
            .filter(|&reaction| {reaction.reaction_type() & 0x1u8 == 0x1u8})
            .collect()
    }
    fn add_element_from_key(&mut self, key:u16, n: usize) {
        for _i in 0..n {
            let form = (key & 0xff00 >> 8) as u8;
            let state = (key & 0xff) as u8;
            let id = self.atoms.len() as u32;
            self.atoms.push(Atom::new(id,form,state));
        }
        self.key_count_cache.clear();
        let x = self.get_atom_number_from_key(key);
    }
    pub fn get_atom_number_from_key(&mut self, key:u16) -> usize {
        if let Some(&count) = self.key_count_cache.get(&key) {
            return count;
        }
        let count = self.atoms.iter().filter(|x| x.get_key() == key).count();
        self.key_count_cache.insert(key, count);
        count
    }
    pub fn get_link_number_from_key(&mut self, key:u32) -> usize {
        if let Some(&count) = self.link_count_cache.get(&key) {
            return count;
        }
        let count = self.links.iter().filter(|&x| self.get_link_key(x) == key).count();
        self.link_count_cache.insert(key, count);
        count
    }
    pub fn get_link_key(&self, link:&Link) -> u32 {
        let a = link.links[0] as usize;
        let b = link.links[1] as usize;
        let a_key = self.atoms[a].get_key();
        let b_key = self.atoms[b].get_key();

        (a_key as u32) << 16 | (b_key as u32)
    }
    pub fn update_all_keys(&mut self) {

    }

    pub fn main_loop(&mut self, tmax: f64, rng: & mut ThreadRng) {
        while self.current_time < tmax {

            let propensity = self.get_propensity_vector();
            let next_time = self.compute_next_reaction_time(&propensity, rng);
            if next_time == None
            {
                // no more reaction
                break;
            }
            else {
                self.current_time += next_time.unwrap();
                let id_reaction = self.compute_next_reaction(&propensity, rng).unwrap();
                let idx = id_reaction;
                //let reaction = & mut self.reactions[idx as usize];
                self.apply_reaction(idx,rng);
            }
            self.cleanup();
        }
    }


    pub fn get_propensity_vector(&mut self) -> Vec<f64> {
        let mut propensity: Vec<f64> = Vec::new();
        for i in 0..self.reactions.len() {
            let reaction:&Reaction = &self.reactions[i];
            if reaction.is_substrate_contact() {
                let substrate_key_a = reaction.get_substrate_a();
                let substrate_key_b = reaction.get_substrate_b();
                // let mut na = self.get_atom_number_from_key(substrate_key_a);
                // let mut nb = self.get_atom_number_from_key(substrate_key_b);
                // let fa = (self.get_atom_number_from_key(substrate_key_a) as f64 / self.nb_current_atoms as f64);
                // let fb = nb as f64 / self.nb_current_atoms as f64;
                propensity[i] = reaction.reaction_speed()
                    * (self.get_atom_number_from_key(substrate_key_a) as f64 / self.nb_current_atoms as f64)
                    * (self.get_atom_number_from_key(substrate_key_b) as f64 / self.nb_current_atoms as f64);
            }
            else {
                let substrate_key_a = reaction.get_substrate_a();
                let substrate_key_b = reaction.get_substrate_b();
                let main_key = (substrate_key_a as u32) << 16 | (substrate_key_b as u32);
                propensity[i] =  reaction.reaction_speed()
                    * self.get_link_number_from_key(main_key) as f64
            }
        }
        propensity
    }

    pub fn compute_a0(& self, propensity : &Vec<f64>) -> f64 {
        let mut a0 = 0f64;

        for i in 0..propensity.len() {
            a0 += propensity[i];
        }
        a0
    }
    pub fn compute_next_reaction_time(&mut self, propensity: &Vec<f64>, rng:&mut ThreadRng) -> Option<f64>{

        let a0 = self.compute_a0(propensity);
        if a0 < 1e-8 {
            return None;
        }
        let mut o = rng.random_range(0f64..1f64);
        Some(-o.ln()/a0)
    }
    pub fn compute_next_reaction(&mut self, propensity: &Vec<f64>, rng:&mut ThreadRng) -> Option<u32> {
        let mut idx = 0i32;
        let a0  = self.compute_a0(propensity);
        let mut o = rng.random_range(0f64..1f64);
        let mut r=0f64;
        while r <= o  && idx < propensity.len() as i32 {
            r += propensity[idx as usize] / a0;
            idx += 1;
        }
        if idx >= propensity.len() as i32 {
            return None;
        }
        idx += -1;
        Some(idx as u32)
        //Some(&self.reactions[idx as usize])
    }

    pub fn apply_reaction(&mut self, idx: u32, rng:&mut ThreadRng) {
        let reaction = &self.reactions[idx as usize] as &Reaction;
        if reaction.is_substrate_contact() {
            // collision type interaction for this reaction
            let substrate_key_a = reaction.get_substrate_a();
            let substrate_key_b = reaction.get_substrate_b();
            if let Some(ida) = self.get_random_atom_id(substrate_key_a, rng) {
                if let Some(idb) = self.get_random_atom_id(substrate_key_b, rng) {
                    let product_state_a = (reaction.get_product_a() & 0x00FF) as u8;
                    let product_state_b = (reaction.get_product_b() & 0x00FF) as u8;
                    self.atoms[ida as usize].set_state(product_state_a);
                    self.atoms[idb as usize].set_state(product_state_b);

                    if !reaction.is_product_contact(){
                        // collision creates a link
                        // example a+b -> c.d
                        self.atoms[ida as usize].link(idb as u32);
                        self.atoms[idb as usize].link(ida as u32);
                        self.links.push(Link{id: self.links.len() as u32,
                            links: [ida as u32, idb as u32],
                            destroy:false});
                    } else {
                        // collision does not create a link
                        // example a+b -> c + d
                        // do nothing
                    }
                } else {
                    panic!("NO REACTANT FOUND shouldn't happen");
                }
            } else {
                panic!("NO REACTANT FOUND shouldn't happen");}
            }
        else {
            // link style reaction
            let substrate_key_a = reaction.get_substrate_a();
            let substrate_key_b = reaction.get_substrate_b();
            let key1 = (substrate_key_a as u32) << 16 | (substrate_key_b as u32);
            let link = self.get_random_link(key1, rng).unwrap();
            if let Some(link) = self.get_random_link(key1, rng) {
                // apply product reaction
                let ida = link.links[0] as usize;
                let idb = link.links[1] as usize;
                let product_state_a = (reaction.get_product_a() & 0x00FF) as u8;
                let product_state_b = (reaction.get_product_b() & 0x00FF) as u8;
                if reaction.is_product_contact() {
                    // collision unlink
                    // example a+b -> c + d
                    self.link_remove.push(link.id);
                    self.atoms[ida as usize].unlink(idb as u32);
                    self.atoms[idb as usize].unlink(ida as u32);

                } else {
                    // collision does not remove a link
                    // example a+b -> c.d
                    // do nothing
                }
                self.atoms[ida as usize].set_state(product_state_a);
                self.atoms[idb as usize].set_state(product_state_b);


            } else {
                    panic!("NO LINK FOUND shouldn't happen");
            }
        }
    }

    pub fn cleanup(&mut self) {

        for idx in 0..self.link_remove.len() {
            self.links.remove(self.link_remove[idx] as usize);
        }
        self.link_remove.clear();
    }

    pub fn get_random_atom_id(&self, key:u16, rng:&mut ThreadRng) -> Option<i32> {
        let v : Vec<&Atom> = self.atoms.iter()
            .filter(|x| x.get_key() == key)
            .collect();
        if v.len() == 0 {
            return None
        }
        else {
            let idx = rng.random_range(0..v.len());
            Some(v[idx].id() as i32)
        }
    }
    pub fn get_random_link(&self, key:u32, rng:&mut ThreadRng) -> Option<&Link> {
        let v : Vec<&Link> = self.links.iter()
            .filter(|x| self.get_link_key(x) == key)
            .collect();
        if v.len() != 0 {
            let idx = rng.random_range(0..v.len());
            Some(v[idx])
        }
        else {
            None
        }
    }
}