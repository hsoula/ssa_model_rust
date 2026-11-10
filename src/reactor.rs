use crate::atoms::Atom;
use crate::links::{get_key, Link};
use crate::event_log::{EventLog};
use crate::molecules::Molecule;
use crate::reactions::{Reaction, ReactionFlag};
use std::collections::{HashMap, HashSet, VecDeque};
use std::process::id;
use rand::prelude::StdRng;
use rand::Rng;
use crate::trace;
use crate::utils::{ObjectId};
use crate::reactions_maps::get_reaction_id;

pub struct Reactor {
    links: HashMap<ObjectId, Link>,
    atoms: HashMap<ObjectId, Atom>,
    reactions: Vec<Reaction>,
    current_time : f64,
    event_log : EventLog,
    nb_link_max: u32,
    fire_counts: Vec<u64>,  // parallel to reactions: how many times each fired
    reaction_sequence: Vec<(f64, usize)>,  // ordered (time, reaction_idx) for every fired reaction
    /// Atom-count snapshot (key → count) and link-count snapshot (key → count)
    /// taken immediately *after* each fired reaction.  Parallel to reaction_sequence.
    state_snapshots: Vec<(HashMap<u16, usize>, HashMap<u32, usize>)>,
}

impl Reactor {
    pub fn new() -> Reactor {
        Reactor {
            links: HashMap::new(),
            atoms: HashMap::new(),
            reactions: Vec::new(),
            current_time: 0f64,
            event_log: EventLog::new(),
            nb_link_max: 4,
            fire_counts: Vec::new(),
            reaction_sequence: Vec::new(),
            state_snapshots: Vec::new(),
        }
    }
    // BOOKMARK: Atoms handling
    // region Atoms
    pub fn add_atom(&mut self, atom: Atom) -> ObjectId {
        let id = atom.id;
        self.atoms.insert(id, atom);
        id
    }
    // Suppression
    pub fn remove_atom(&mut self, id: ObjectId)  {
        let atom = self.atoms.get_mut(&id).unwrap();
        let mut linked = atom.unlink_all();
        for link in linked.iter() {
            self.get_atom_mut(*link).unwrap().unlink(id);
        }
        let remove : Vec<ObjectId> = self.links.iter().filter(|(idx, l)| l.links[0]==id || l.links[1]==id)
            .map(|(idx, l)| *idx)
            .collect();
        for r in remove { self.links.remove(&r); }
        self.atoms.remove(&id);
    }

    // Find
    pub fn get_atom(&self, id: ObjectId) -> Option<&Atom> {
        self.atoms.get(&id)
    }

    // mutable find
    pub fn get_atom_mut(&mut self, id: ObjectId) -> Option<&mut Atom> {
        self.atoms.get_mut(&id)
    }

    // contains - O(1)
    pub fn contains_atom(&self, id: ObjectId) -> bool {
        self.atoms.contains_key(&id)
    }

    //
    pub fn get_atom_number(&self) -> usize {
        self.atoms.len()
    }

    pub fn get_atom_with_max_links(&self) -> usize {
        self.atoms.iter().filter(|(id,atom)| atom.get_nb_link() >= self.nb_link_max ).count()
    }
    pub fn get_atom_number_from_key(&self, key:u16) -> usize {
        self.get_atom_from_key(key).len()
    }

    // Iter over all atoms
    pub fn iter(&self) -> impl Iterator<Item = &Atom> {
        self.atoms.values()
    }

    // Itération mutable
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Atom> {
        self.atoms.values_mut()
    }

    pub fn get_atom_from_key(&self, key:u16) -> Vec<&Atom> {
        self.atoms.values()
           // .filter(|atom| ((atom.get_key() & key ) & 0xff00 == atom.get_key() & 0xff00 ) && (atom.get_key() & 0x00ff == key & 0x00ff))
            .filter(|atom| atom.get_key() == key)
            .collect()
    }
    //  clear atoms
    pub fn clear_atoms(&mut self) {
        self.atoms.clear();
    }

    // vector of ids
    pub fn get_all_ids(&self) -> Vec<ObjectId> {
        self.atoms.keys().cloned().collect()
    }
    pub fn add_atom_from_keys(&mut self, key:u16, n: usize) {
        for _i in 0..n {

            self.add_atom(Atom::new_from_key(key));
        }
    }
    pub fn get_random_atom_id(&self, key:u16, rng:&mut StdRng) -> Option<ObjectId> {
        let v : Vec<&Atom> = self.get_atom_from_key(key);
        if v.len() == 0 {
            return None
        }
        else {
            let idx = rng.random_range(0..v.len());
            Some(v[idx].id() as ObjectId)
        }
    }
    // endregion

    // BOOKMARK: reaction handling
    /// Reaction functions
    ///

    pub fn add_reaction(&mut self, reaction: Reaction) {
        self.reactions.push(reaction);
        self.fire_counts.push(0);
    }

    pub fn get_reactions_key(&self, idx:usize) -> u32 {
        self.reactions[idx].substrate_key()
    }

    pub fn get_collision_reactions(&self) -> Vec<&Reaction> {
        self.reactions.iter()
            .filter(|&reaction| {reaction.reaction_type() & (ReactionFlag::ProductContact as u8) == (ReactionFlag::ProductContact as u8)})
            .collect()
    }


    /// links function
    ///
    pub fn add_link(&mut self, link: Link) -> ObjectId {
        let id = link.id;
        self.links.insert(id, link);
        id
    }
    // Suppression
    pub fn remove_link(&mut self, id: ObjectId) -> Option<Link> {
        self.links.remove(&id)
    }

    // Find
    pub fn get_link(&self, id: ObjectId) -> Option<&Link> {
        self.links.get(&id)
    }

    // mutable find
    pub fn get_link_mut(&mut self, id: ObjectId) -> Option<&mut Link> {
        self.links.get_mut(&id)
    }

    // contains - O(1)
    pub fn contains_link(&self, id: ObjectId) -> bool {
        self.links.contains_key(&id)
    }

    //
    pub fn get_link_number(&self) -> usize {
        self.links.len()
    }

    pub fn get_link_number_from_key(&self, key:u32) -> usize {
        self.get_link_from_key(key).len()
    }

    // Iter over all atoms
    pub fn iter_links(&self) -> impl Iterator<Item = &Link> {
        self.links.values()
    }

    // Itération mutable
    pub fn iter_mut_links(&mut self) -> impl Iterator<Item = &mut Link> {
        self.links.values_mut()
    }

    pub fn get_link_from_key(&self, key:u32) -> Vec<&Link> {
        self.links.values()
            .filter(|link| link.key() == key || link.key_rev() == key)
            .collect()
    }
    //  clear atoms
    pub fn clear_links(&mut self) {
        self.links.clear();
    }

    // vector of ids
    pub fn get_all_links_ids(&self) -> Vec<ObjectId> {
        self.links.keys().cloned().collect()
    }
    pub fn get_random_link_id(&self, key:u32, rng:&mut StdRng) -> Option<ObjectId> {
        let v : Vec<&Link> = self.get_link_from_key(key);
        if v.len() == 0 {
            return None
        }
        else {
            let idx = rng.random_range(0..v.len());
            Some(v[idx].id() as ObjectId)
        }
    }

    pub fn main_loop(&mut self, tmax: f64, rng: & mut StdRng) -> f64{
        while self.current_time < tmax {
            // compute propensity
            let propensity = self.get_propensity_vector();
            // random select the next reaction time
            let next_time = self.compute_next_reaction_time(&propensity, rng);

            // no reactions possible
            if next_time == None
            {
                // no more reaction
                break;
            } else {
                self.current_time += next_time.unwrap();
                // select the reaction
                let id_reaction = self.compute_next_reaction(&propensity, rng).unwrap();
                let idx = id_reaction;
                let key = self.get_reactions_key(idx as usize);

                // apply the reaction e.g modify the molecule graphs
                self.apply_reaction(idx, rng);
                if let Some(c) = self.fire_counts.get_mut(idx as usize) { *c += 1; }
                self.reaction_sequence.push((self.current_time, idx as usize));
                // Snapshot: atom-counts and link-counts right after this reaction.
                let atom_snap: HashMap<u16, usize> = self.atoms.values()
                    .map(|a| a.get_key())
                    .fold(HashMap::new(), |mut acc, k| { *acc.entry(k).or_insert(0) += 1; acc });
                let link_snap: HashMap<u32, usize> = self.links.values()
                    .map(|l| l.key())
                    .fold(HashMap::new(), |mut acc, k| { *acc.entry(k).or_insert(0) += 1; acc });
                self.state_snapshots.push((atom_snap, link_snap));

                // save reaction events in cas they have tags
                if let Some(tag) = self.reactions[id_reaction as usize].tags() {
                    if let Some(u) = get_reaction_id(tag.as_str()) {
                        self.event_log.add_event(self.current_time, u);
                    }
                }
            }
            // self.display_stats();
            // self.cleanup();
        }
        self.current_time
    }

    pub fn compute_mask_reaction_number(&self, sub1_key:u16, sub2_key:u16) -> f64 {
        let mut a = 0f64;
        let sub1_mask = ((sub1_key & 0xff00) >> 8) as u8;
        let sub2_mask = ((sub2_key & 0xff00) >> 8) as u8;
        if sub1_mask == sub2_mask {

            // difficult mode  we affect the all product from every pair
            let sub1_stats = self.atoms.values()
                .map(|atom| atom.get_key())
                .filter(|&atom_key| (atom_key & sub1_key) == atom_key)
                .fold(HashMap::new(), |mut acc, key| {
                    *acc.entry(key).or_insert(0) += 1;
                    acc
                });
            let x = self.atoms.values()
                .filter(|&atom| (atom.get_key() & sub1_key) == atom.get_key()).count() ;
            if sub1_key == sub2_key {

                for (_key, value) in sub1_stats.iter() {
                    let v = *value as f64;
                    a += v * (v - 1f64);
                }
            } else {

                let sub2_stats = self.atoms.values()
                    .map(|atom| atom.get_key())
                    .filter(|&atom_key| (atom_key & sub2_key) == atom_key)
                    .fold(HashMap::new(), |mut acc, key| {
                        *acc.entry(key).or_insert(0) += 1;
                        acc
                    });

                for (key, value) in sub1_stats.iter() {
                    let key_mask = ((key & 0xff00) >> 8) as u8;
                    let v = (*value) as f64;
                    let c :HashMap<&u16,&i32> = sub2_stats.iter().
                        filter(|&(key, _)| ((key & 0xff00) >> 8) as u8 == key_mask)
                        .collect();
                    let x = c.iter()
                        .map(|(&key, &value)| *value as f64).sum::<f64>();
                    a += v * x;
                }
            }
        }
        else {
                // sub1_mask != sub2_mask
                let a1 = self.get_atom_number_from_key(sub1_key) as f64;
                let b1 = self.get_atom_number_from_key(sub2_key) as f64;

                a = a1 * b1;
            }
        a
    }

    pub fn get_propensity_vector(&mut self) -> Vec<f64> {
        let mut propensity: Vec<f64> = vec![0.0; self.reactions.len()];
        for i in 0..self.reactions.len() {
            let reaction:&Reaction = &self.reactions[i];

            if reaction.is_production() {
                propensity[i] = reaction.reaction_speed(); // lol easy ^^
            }
            else if reaction.is_degradation() {
                let substrate_key_a = reaction.get_substrate_a();
                propensity[i] = reaction.reaction_speed() * (self.get_atom_number_from_key(substrate_key_a) as f64);
            }
            else if reaction.is_substrate_contact() {
                let substrate_key_a = reaction.get_substrate_a();
                let substrate_key_b = reaction.get_substrate_b();

                // let mut na = self.get_atom_number_from_key(substrate_key_a);
                // let mut nb = self.get_atom_number_from_key(substrate_key_b);
                // let fa = (self.get_atom_number_from_key(substrate_key_a) as f64 / self.nb_current_atoms as f64);
                // let fb = nb as f64 / self.nb_current_atoms as f64;
                let mut offset = 0f64;
                if substrate_key_a == substrate_key_b {
                    offset = -1f64;
                }
                let key_link = (substrate_key_a as u32) << 16 | (substrate_key_b as u32);
                let linked = self.get_link_from_key(key_link).len() as f64 ;
                let key_link_rev = (substrate_key_b as u32) << 16 | (substrate_key_a as u32);
                let linked_rev = self.get_link_from_key(key_link_rev).len() as f64 ;
                let key_a = self.get_atom_number_from_key(substrate_key_a);
                let key_b = self.get_atom_number_from_key(substrate_key_b);
                let s = reaction.to_string();
                propensity[i] = reaction.reaction_speed() //* self.compute_mask_reaction_number(substrate_key_a, substrate_key_b);
                    * ((self.get_atom_number_from_key(substrate_key_a) as f64 - linked - linked_rev).max(0f64) //    / self.get_atom_number() as f64)
                    * (self.get_atom_number_from_key(substrate_key_b) as f64 + offset - linked_rev - linked).max(0f64)); // / self.get_atom_number() as f64);
            }
            else {
                let substrate_key_a = reaction.get_substrate_a();
                let substrate_key_b = reaction.get_substrate_b();

                let main_key = (substrate_key_a as u32) << 16 | (substrate_key_b as u32);
                let main_key_rev = (substrate_key_b as u32) << 16 | (substrate_key_a as u32);
                propensity[i] =  reaction.reaction_speed()
                    * (self.get_link_number_from_key(main_key) as f64 + self.get_link_number_from_key(main_key_rev) as f64);
            }
        }
     //   for l in self.links.values() {
      //  }

        propensity
    }

    pub fn compute_a0(& self, propensity : &Vec<f64>) -> f64 {
        let mut a0 = 0f64;

        for i in 0..propensity.len() {
            a0 += propensity[i];
        }
        a0
    }
    pub fn compute_next_reaction_time(&mut self, propensity: &Vec<f64>, rng:&mut StdRng) -> Option<f64>{

        let a0 = self.compute_a0(propensity);
        if a0 < 1e-8 { // TODO put a general simulation constant
            return None;
        }
        let mut o = rng.random_range(0f64..1f64);
        Some(-o.ln()/a0)
    }
    pub fn compute_next_reaction(&mut self, propensity: &Vec<f64>, rng:&mut StdRng) -> Option<u32> {
        let mut idx = 0i32;
        let a0  = self.compute_a0(propensity);
        let o = rng.random_range(0f64..1f64);
        let mut r=0f64;
        while r <= o  && idx < propensity.len() as i32 {
            r += propensity[idx as usize] / a0;
            idx += 1;
        }
        if idx >= propensity.len() as i32 + 1 {
            return None;
        }
        idx += -1;
        Some(idx as u32)
    }

    pub fn apply_contact_reaction(&mut self, ida : ObjectId, idb : ObjectId, product_state_a: u8, product_state_b: u8, product_contact: bool) {
        self.get_atom_mut(ida).unwrap().set_state(product_state_a);
        self.get_atom_mut(idb).unwrap().set_state(product_state_b);


        if !product_contact{
            // collision creates a link
            // example a+b -> c.d
            if  self.get_atom_mut(ida).unwrap().get_nb_link() < self.nb_link_max
                && self.get_atom_mut(idb).unwrap().get_nb_link() < self.nb_link_max {
                self.get_atom_mut(ida).unwrap().link(idb);
                self.get_atom_mut(idb).unwrap().link(ida);
                let mut link = Link::new(ida, self.get_atom(ida).unwrap().get_key(), idb, self.get_atom(idb).unwrap().get_key());
                self.add_link(link);
            }
        } else {
            // collision does not create a link
            // example a+b -> c + d
            // do nothing
        }
    }

    pub fn apply_link_reaction(&mut self, link_id : ObjectId, product_state_a : u8, product_state_b:u8, product_contact: bool) {
        // apply product reaction
        let link = self.get_link(link_id).unwrap();
        let ida = link.links[0];
        let idb = link.links[1];
        if product_contact {
            // collision unlink
            // example a+b -> c + d
            self.remove_link(link.id);
            self.get_atom_mut(ida).unwrap().unlink(idb);
            self.get_atom_mut(idb).unwrap().unlink(ida);

        } else {
            // collision does not remove a link
            // example a+b -> c.d
            // do nothing
        }

        self.get_atom_mut(ida).unwrap().set_state(product_state_a);
        self.get_atom_mut(idb).unwrap().set_state(product_state_b);
    }

    pub fn get_random_atoms_from_reaction (&self, substrate_key_a:u16, substrate_key_b:u16,rng: &mut StdRng) -> (ObjectId, ObjectId) {
        let sub1_mask = ((substrate_key_a & 0xff00) >> 8) as u8;
        let sub2_mask = ((substrate_key_b & 0xff00) >> 8) as u8;
        let ida = self.get_random_atom_id(substrate_key_a, rng).unwrap();
        if sub1_mask == sub2_mask {
            // difficult mode  we affect the all product from every pair
            // we need to find an atom from the same form
            // TODO we need to check if there is a possible pair
            let new_mask = (self.get_atom(ida).unwrap().get_key() & 0xFF00) | (substrate_key_b & 0x00FF);
            let idb = self.get_random_atom_id(new_mask, rng).unwrap();

            (ida,idb)
            }
        else {
            let idb = self.get_random_atom_id(substrate_key_b, rng).unwrap();

            (ida, idb)
        }
    }
    pub fn apply_reaction(&mut self, idx: u32, rng:&mut StdRng) {
        let reaction = &self.reactions[idx as usize] as &Reaction;
        let substrate_key_a = reaction.get_substrate_a();
        let substrate_key_b = reaction.get_substrate_b();
        let product_state_a = (reaction.get_product_a() & 0x00FF) as u8;
        let product_state_b = (reaction.get_product_b() & 0x00FF) as u8;
        let substrate_contact = reaction.is_substrate_contact();
        let product_contact = reaction.is_product_contact();
        let s = reaction.to_string();
        let key_a =substrate_key_a;
        let key_b =substrate_key_b;
        if reaction.is_degradation() {
            let id = self.get_random_atom_id(substrate_key_a, rng).unwrap();
            self.remove_atom(id);
        } else if reaction.is_production() {
            // add one atome for the reaction
            self.add_atom_from_keys(substrate_key_a, 1);
        } else {

            if substrate_contact {
                // collision type interaction for this reaction
                // we need to find guys not linked
                let mut done = false;
                while !done {
                    let (ida, idb) = self.get_random_atoms_from_reaction(substrate_key_a, substrate_key_b, rng);
                    if (ida == idb) {continue;}
                    if !self.atoms.get(&ida).unwrap().is_link(idb) {
                        self.apply_contact_reaction(ida, idb, product_state_a, product_state_b, product_contact);
                        done = true;
                    }
                }
            } else {
                // link style reaction
                let key1 = (substrate_key_a as u32) << 16 | (substrate_key_b as u32);
                if let Some(link_id) = self.get_random_link_id(key1, rng) {
                    self.apply_link_reaction(link_id, product_state_a, product_state_b, product_contact);
                }
            }
        }
        self.relink();

    }

    pub fn relink(&mut self) {

        let updates: Vec<(ObjectId, u16, u16)> = self.links.iter()
            .map(|(&link_id, link)| {
                let atom_a = self.get_atom(link.links[0]).unwrap();
                let atom_b = self.get_atom(link.links[1]).unwrap();
                (link_id, atom_a.get_key(), atom_b.get_key())
            })
            .collect();

        // Puis appliquer les modifications
        for (link_id, ka, kb) in updates {
            if let Some(link) = self.links.get_mut(&link_id) {
                link.update_key(ka, kb);
            }
        }
    }

    pub fn display_stats(&self) -> (HashMap<String,i32>, HashMap<u32, i32>) {
        let state_states = self.atoms.values()
            .map(|atom| atom.to_string())
            .fold(HashMap::new(), |mut acc, key| {
                *acc.entry(key).or_insert(0) += 1;
                acc
            });
        let link_states = self.links.values()
            .map(|link| link.key())
            .fold(HashMap::new(), |mut acc, key| {
                *acc.entry(key).or_insert(0) += 1;
                acc
            });

        // let form_states = self.atoms.values()
        //     .map(|atom| atom.form())
        //     .fold(HashMap::new(), |mut acc, key| {
        //         *acc.entry(key).or_insert(0) += 1;
        //         acc
        //     });

        for (key, value) in state_states.iter() {
            println!("{} -> {}", key, value);
        }
        for (key, value) in link_states.iter() {
            println!("{:08x} -> {}", key, value);
        }
        (state_states, link_states)
    }

    pub fn display_events(&self) -> HashMap<String,i32> {
       let event = self.event_log.event_stats();
        for (key, value) in event.iter() {
            println!("{:} -> {:}", key, value);
        }
        event
    }

    pub fn add_molecule(&mut self, molecule: &Molecule) {

        let mut index = HashMap::new();
        for (key, value) in molecule.atoms.iter() {
            let form = value.0;
            let state = value.1;
            let mut a = Atom::new(form,state);
            let id = a.id();
            self.add_atom(a);
            index.insert(key, id);
        }

        for (key, value) in molecule.links.iter() {
            let key1 = value.0;
            let key2 = value.1;
            let o1 = index.get(&key1).unwrap().clone();
            let o2 = index.get(&key2).unwrap().clone();
            self.get_atom_mut(o1).unwrap().link(o2);
            self.get_atom_mut(o2).unwrap().link(o1);
            let mut link = Link::new(o1,self.get_atom(o1).unwrap().get_key(), o2,self.get_atom(o2).unwrap().get_key());
            self.add_link(link);
        }

    }

    pub fn export_molecule(&self) -> Molecule {
        let mut molecule = Molecule::new();
        let mut index = HashMap::new();

        for (key, value) in self.atoms.iter().enumerate() {
            index.insert(value.0 , key as u32);
            molecule.atoms.insert(key as u32, (value.1.form(), value.1.state()));
        }
        for (key, value) in self.links.iter().enumerate() {
            let key1 = &value.1.links[0];
            let idx1 = index.get(key1).unwrap();
            let key2 = &value.1.links[1];
            let idx2 = index.get(key2).unwrap();

            molecule.links.insert(key as u32 + 1,(idx1.clone(), idx2.clone()));
        }
        molecule
    }

    pub fn export_to_file(&self, file: &str) {
        let molecule = self.export_molecule();
        molecule.save_to_file(file);
    }


    pub fn get_events_number(&self) {}

    pub fn reactions(&self) -> &[Reaction] { &self.reactions }
    pub fn current_time(&self) -> f64 { self.current_time }
    pub fn reaction_sequences(&self) -> &[(f64, usize)] { &self.reaction_sequence }
    /// Per-step snapshots: atom-count map and link-count map, parallel to `reaction_sequences()`.
    pub fn state_snapshots(&self) -> &[(HashMap<u16, usize>, HashMap<u32, usize>)] {
        &self.state_snapshots
    }

    /// Returns (reaction, fire_count) pairs sorted by fire count descending.
    pub fn reaction_fire_counts(&self) -> Vec<(&Reaction, u64)> {
        let mut v: Vec<(&Reaction, u64)> = self.reactions.iter()
            .zip(self.fire_counts.iter().copied())
            .collect();
        v.sort_by(|a, b| b.1.cmp(&a.1));
        v
    }

    // ── Markdown report ──────────────────────────────────────────────────────────

    /// Render the full reactor state as a Markdown string.
    pub fn to_markdown(&self) -> String {
        let mut md = String::new();

        md.push_str(&format!("# Reactor State — t = {:.3}\n\n", self.current_time));

        // ── Atoms ──
        let atom_counts: HashMap<String, usize> = self.atoms.values()
            .map(|a| a.to_string())
            .fold(HashMap::new(), |mut acc, k| { *acc.entry(k).or_insert(0) += 1; acc });

        md.push_str(&format!("## Atoms ({} total)\n\n", self.atoms.len()));
        md.push_str("| Symbol | Count |\n|--------|-------|\n");
        let mut asym: Vec<_> = atom_counts.iter().collect();
        asym.sort_by_key(|(k, _)| k.clone());
        for (sym, count) in &asym {
            md.push_str(&format!("| {:<6} | {:>5} |\n", sym, count));
        }

        // ── Links ──
        let link_counts: HashMap<u32, usize> = self.links.values()
            .map(|l| l.key())
            .fold(HashMap::new(), |mut acc, k| { *acc.entry(k).or_insert(0) += 1; acc });

        md.push_str(&format!("\n## Links ({} total)\n\n", self.links.len()));
        if link_counts.is_empty() {
            md.push_str("*(none)*\n");
        } else {
            md.push_str("| Bond        | Count |\n|-------------|-------|\n");
            let mut lsym: Vec<_> = link_counts.iter().collect();
            lsym.sort_by_key(|(k, _)| *k);
            for (key, count) in &lsym {
                let a = reactor_key_sym((**key >> 16) as u16);
                let b = reactor_key_sym((**key & 0xffff) as u16);
                md.push_str(&format!("| {}—{} | {:>5} |\n", a, b, count));
            }
        }

        // ── Reactions ──
        md.push_str(&format!("\n## Reactions ({} total)\n\n", self.reactions.len()));
        for r in &self.reactions {
            let rtype = reaction_type_label(r);
            md.push_str(&format!(
                "- `{}` *({}, speed={:.3})*\n",
                r.to_string(), rtype, r.reaction_speed()
            ));
        }

        // ── Molecule size distribution ──
        let mols = self.find_molecules();
        let mut size_dist: HashMap<usize, usize> = HashMap::new();
        for m in &mols { *size_dist.entry(m.size()).or_insert(0) += 1; }

        md.push_str(&format!("\n## Molecules ({} connected components)\n\n", mols.len()));
        md.push_str("| Size | Count |\n|------|-------|\n");
        let mut sizes: Vec<_> = size_dist.iter().collect();
        sizes.sort_by_key(|(k, _)| *k);
        for (size, count) in &sizes {
            md.push_str(&format!("| {:>4} | {:>5} |\n", size, count));
        }

        // ── Reaction fire counts ──
        let total_fires: u64 = self.fire_counts.iter().sum();
        md.push_str(&format!("\n## Reaction fire counts ({} total)\n\n", total_fires));
        md.push_str("| Rank | Reaction | Fires | % |\n|------|----------|-------|---|\n");
        let by_fire = self.reaction_fire_counts();
        for (rank, (r, count)) in by_fire.iter().enumerate() {
            let pct = if total_fires > 0 { 100.0 * (*count) as f64 / total_fires as f64 } else { 0.0 };
            md.push_str(&format!("| {:>4} | `{}` | {:>6} | {:>5.1}% |\n",
                rank + 1, r.to_string(), count, pct));
        }

        md
    }

    /// Write the markdown report to a file.
    pub fn save_markdown(&self, path: &str) {
        std::fs::write(path, self.to_markdown())
            .unwrap_or_else(|e| eprintln!("Cannot write markdown to {path}: {e}"));
    }

    pub fn find_molecules(&self) -> Vec<Molecule> {

        let mut visited = HashSet::new();
        let mut molecules = Vec::new();

        for atom in self.iter() {
            if !visited.contains(&atom.id) {

                let mut molecule = Molecule::new();
                let mut queue = VecDeque::new();
                let mut component_atoms = HashSet::new();

                queue.push_back(atom.id);
                visited.insert(atom.id);
                component_atoms.insert(atom.id);

                while let Some(current_id) = queue.pop_front() {
                    let current_atom = self.get_atom(current_id).unwrap();
                    let state = current_atom.state();
                    let form = current_atom.form();

                    molecule.add_atom(current_id as u32, (form, state));


                    for &neighbor_id in current_atom.get_links() {
                        let key1 = self.get_atom(current_id).unwrap().get_key();
                        let key2 = self.get_atom(neighbor_id).unwrap().get_key();
                        //let link = Link::new(current_id, key1, neighbor_id, key2);

                        if !molecule.has_link((current_id as u32, neighbor_id as u32)) && !molecule.has_link((neighbor_id as u32, current_id as u32)){
                            let id = molecule.links.len() as u32;
                            molecule.links.insert(id, (current_id as u32, neighbor_id as u32));
                        }

                        // Ajouter les voisins non visités
                        if !visited.contains(&neighbor_id) {
                            visited.insert(neighbor_id);
                            component_atoms.insert(neighbor_id);
                            queue.push_back(neighbor_id);
                        }
                    }
                }

                molecules.push(molecule);
            }
        }

        molecules
    }
}

/// Convert an atom u16 key (form<<8 | state) to its symbolic label, e.g. "a0", "b3".
fn reactor_key_sym(key: u16) -> String {
    let form  = ((key & 0xff00) >> 8) as u8;
    let state = (key & 0x00ff) as u8;
    format!("{}{}", (form + 96) as char, state)
}

fn reaction_type_label(r: &Reaction) -> &'static str {
    if r.is_production()  { return "production"; }
    if r.is_degradation() { return "degradation"; }
    match (r.is_substrate_contact(), r.is_product_contact()) {
        (true,  true)  => "contact→bond",
        (true,  false) => "contact",
        (false, true)  => "→bond",
        (false, false) => "state-change",
    }
}