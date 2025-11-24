use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use serde_with::{serde_as, DisplayFromStr};
use crate::utils::JsonReaction;
//#[derive(Debug, Clone, Copy, Deserialize)]

#[repr(u8)]
pub enum ReactionFlag {
    SubstrateContact = 0b10u8,
    ProductContact = 0b1u8,
}
//pub mod flags;
#[derive(Debug)]
pub struct Reaction {
    id: u32,
    reaction_type: u8,
    reaction_speed : f64,
    substrate_key: u32,
    product_key: u32
}

impl Reaction {
    pub fn new(id: u32) -> Reaction {
        Reaction{id:id,
            reaction_type: 0,
            reaction_speed: 1.0,
            substrate_key:0,
            product_key:0}
    }

    pub fn new_from_json(json: &str, id: u32) -> Reaction {
        let file = File::open(json).unwrap();
        let reader = BufReader::new(file);
        let js :JsonReaction =  serde_json::from_reader(reader).unwrap();
        let rt = js.reaction_type as u8;
        let rs = js.reaction_speed as f64;
        let sk = (js.s1 as u32) << 16 | (js.s2 as u32);
        let pk = (js.p1 as u32) << 16 | (js.p2 as u32);
        Reaction{id:id, reaction_type:rt, reaction_speed:rs, substrate_key:sk, product_key: pk}

    }
    pub fn reaction_type(&self) -> u8 {
        self.reaction_type
    }

    pub fn is_substrate_contact(&self) -> bool {
        self.reaction_type & (ReactionFlag::SubstrateContact as u8) == (ReactionFlag::SubstrateContact as u8)
    }
    pub fn is_product_contact(&self) -> bool {
        self.reaction_type & (ReactionFlag::ProductContact as u8) == (ReactionFlag::ProductContact as u8)
    }

    pub fn reaction_speed(&self) -> f64 {
        self.reaction_speed
    }

    pub fn get_substrate_a(&self) -> u16 {
        ((self.substrate_key & 0xffff0000) >> 16) as u16
    }
    pub fn get_substrate_b(&self) -> u16 {
        (self.substrate_key & 0x0000ffff)  as u16
    }
    pub fn get_product_a(&self) -> u16 {
        ((self.product_key & 0xffff0000) >> 16) as u16
    }
    pub fn get_product_b(&self) -> u16 {
        (self.product_key & 0x0000ffff)  as u16
    }
}