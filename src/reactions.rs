//use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use std::option::Option;
//use serde::de::Unexpected::Option;
use crate::trace;
//use serde_with::{serde_as, DisplayFromStr};
use crate::utils::JsonReaction;
use crate::atoms::{key_symbol, OUTSIDE};
//#[derive(Debug, Clone, Copy, Deserialize)]
//use crate::trace;

#[repr(u8)]
pub enum ReactionFlag {
    SubstrateContact = 0b10u8,
    ProductContact = 0b01u8,
    Production = 0b100u8,
    Degradation = 0b1000u8,
    /// Unimolecular state change of a single atom: `a1 -> a11`.
    /// Used e.g. for membrane transport in the bacterium reactor.
    Transform = 0b10000u8,
}
//pub mod flags;
#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    id: u32,
    reaction_type: u8,
    reaction_speed : f64,
    substrate_key: u32,
    product_key: u32,
    tags: Option<String>
}

impl Reaction {
    pub fn new(id: u32) -> Reaction {
        Reaction{id:id,
            reaction_type: 0,
            reaction_speed: 1.0,
            substrate_key:0,
            product_key:0,
            tags:None,}
    }
    /*
        start : 0
        end : 1
        sf1 : 2
        ss1 : 3
        sf2 : 4
        ss2 : 5
        c : 6
        pf1 : 7
        ps1 : 8
        pf2 : 9
        ps2 : 10
        speed : 11
     */
    pub fn new_from_vec8(id: u32, vec: Vec<u8>, max_count: u8) -> Reaction {
        let mut ps1 = 0u8;
        let mut ps2 = 0u8;
        let mut sf1 = 0u8;
        let mut sf2 = 0u8;
        let mut ss1 = 0u8;
        let mut ss2 = 0u8;

        let mut contact = 0u8;
        let mut reaction_speed = 0.1f64; // minimum speed is 0.1; each D token adds 0.1
        let max_form  = max_count;
        let max_state = max_count.saturating_mul(2); // states have 2× the range of forms

        for &token in vec.iter() {
            match token {
                0x0E             => {}                                                   // E = start marker
                0x0F             => {}                                                   // F = end marker
                0x00 | 0x01      => { if sf1 < max_form  { sf1 += 1 } },                // substrate 1 form  (x2)
                0x02 | 0x03      => { if ss1 < max_state { ss1 += 1 } },                // substrate 1 state (x2)
                0x04 | 0x05      => { if sf2 < max_form  { sf2 += 1 } },                // substrate 2 form  (x2)
                0x06 | 0x07      => { if ss2 < max_state { ss2 += 1 } },                // substrate 2 state (x2)
                0x08             => { contact += 1 },                                    // reaction-type bits (uncapped, mod 4 applied below)
                0x09 | 0x0A      => { if ps1 < max_state { ps1 += 1 } },                // product 1 state   (x2)
                0x0B | 0x0C      => { if ps2 < max_state { ps2 += 1 } },                // product 2 state   (x2)
                0x0D             => { reaction_speed += 0.1f64 },                        // speed (additive)
                _                => {}
            }
        }
        contact = contact % 4;
        // Product forms mirror substrate forms — only state can change
        let p1 = ((sf1 as u16) << 8) | (ps1 as u16);
        let p2 = ((sf2 as u16) << 8) | (ps2 as u16);
        let s1 = ((sf1 as u16) << 8) | (ss1 as u16);
        let s2 = ((sf2 as u16) << 8) | (ss2 as u16);
        let substrate_key = ((s1 as u32) << 16) | (s2 as u32);
        let product_key = ((p1 as u32) << 16) | (p2 as u32);
        Reaction{id:id,
        reaction_type:contact,
        reaction_speed:reaction_speed,
        substrate_key: substrate_key,
        product_key:product_key,
        tags:Some("".to_string())}
    }

    pub fn reaction_from_wildcard_all(js : &JsonReaction, offset:u32) -> Vec<Reaction> {
        let mut v = Vec::new();
        let mut ix= offset;
        for i in 0..0x6 {
            for j in 0..0x10 {
                let rt = js.reaction_type;
                let rs = js.reaction_speed;
                let js1 = (i << 8)  | (js.s1 & 0x00FF) as u16;
                let js2 = (j << 8)  | (js.s2 & 0x00FF) as u16;;
                let jp1 = (i << 8)  | (js.p1 & 0x00FF) as u16;
                let jp2 = (j << 8)  | (js.p2 & 0x00FF) as u16;;

                let mut sk = (js1 as u32) << 16 | (js2 as u32);
                let pk = (jp1 as u32) << 16 | (jp2 as u32);
                let s = js.tags.clone();
                v.push( Reaction{ id: ix, reaction_type: rt, reaction_speed: rs, substrate_key: sk, product_key: pk, tags: s});
                ix += 1;
            }
        }
        v
    }
    pub fn reaction_from_wildcard_with(js : &JsonReaction, offset:u32) -> Vec<Reaction> {
        let mut v = Vec::new();
        let mut ix= offset;
        for i in 0..0x10 {
                let rt = js.reaction_type;
                let rs = js.reaction_speed;
                let js1 = (i << 8)  | (js.s1 & 0x00FF) as u16;
                let js2 = (i << 8)  | (js.s2 & 0x00FF) as u16;;
                let jp1 = (i << 8)  | (js.p1 & 0x00FF) as u16;
                let jp2 = (i << 8)  | (js.p2 & 0x00FF) as u16;;


                let mut sk = (js1 as u32) << 16 | (js2 as u32);
                let pk = (jp1 as u32) << 16 | (jp2 as u32);
                let s = js.tags.clone();
                v.push( Reaction{ id: ix, reaction_type: rt, reaction_speed: rs, substrate_key: sk, product_key: pk, tags: s});
                ix += 1;
        }
        v
    }

    pub fn news_from_json(json: &str, id_start: u32) -> Vec<Reaction> {
       // trace!(json);
        let mut v = Vec::new();
        let file = File::open(json).unwrap();
        let reader = BufReader::new(file);
        let jsv :Vec<JsonReaction> =  serde_json::from_reader(reader).unwrap();
        let mut offset = 0;
        for (ix,js) in jsv.iter().enumerate() {
            if js.s1 & 0xFF00 == 0xFF00 && js.s2 & 0xFF00 == 0xFF00 {
                let s = Self::reaction_from_wildcard_all(js, ix as u32);
                offset += s.len();
                v.extend(s);
            }
            else if js.s1 & 0xFF00 == 0xFF00 && js.s2 & 0xEF00 == 0xEF00 {
                let s = Self::reaction_from_wildcard_with(js, ix as u32);
                offset += s.len();
                v.extend(s);
            }
            else {
                let rt = js.reaction_type;
                let rs = js.reaction_speed;
                let sk = (js.s1 as u32) << 16 | (js.s2 as u32);
                let pk = (js.p1 as u32) << 16 | (js.p2 as u32);
                let mut idx = ix as u32 + id_start;
                let s = js.tags.clone();
                v.push(Reaction { id: idx, reaction_type: rt, reaction_speed: rs, substrate_key: sk, product_key: pk, tags: s});
            }
        }
        v
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

    pub fn is_production(&self) -> bool {
        self.reaction_type == ReactionFlag::Production as u8
    }
    pub fn is_degradation(&self) -> bool {
        self.reaction_type == ReactionFlag::Degradation as u8
    }
    pub fn is_transform(&self) -> bool {
        self.reaction_type == ReactionFlag::Transform as u8
    }

    /// Parse a reaction written in the same notation `to_string()` prints.
    ///
    /// Atoms are `<form letter><state>` (a = form 1, b = form 2, ...).
    /// Between two atoms, `+` means free (collision / no bond) and `.` means bonded.
    ///
    /// * `0 -> a1`                 production
    /// * `a1 -> 0`                 degradation
    /// * `a1 -> a11`               transform (single atom changes state)
    /// * `a1 + b12 -> a12.b13`     collision that creates a bond
    /// * `a12.b13 -> a12 + b12`    bonded pair that breaks
    ///
    /// Uppercase letters are outside the bacterium (`A1`), lowercase inside (`a1`).
    ///
    /// Forms never change (Hutton): product i has the same letter as substrate i, only its
    /// case may change (crossing the membrane). Bonded atoms must share a compartment.
    pub fn from_rule(id: u32, rule: &str, speed: f64) -> Result<Reaction, String> {
        let rule: String = rule.chars().filter(|c| !c.is_whitespace()).collect();
        let (lhs, rhs) = rule.split_once("->")
            .ok_or_else(|| format!("missing '->' in rule '{rule}'"))?;
        let base = Reaction::new(id).set_reaction_speed(speed);

        if lhs == "0" {
            let p = parse_atom_token(rhs)?;
            return Ok(base.set_reaction_type(ReactionFlag::Production as u8)
                .set_substrate_key((p as u32) << 16));
        }
        if rhs == "0" {
            let s = parse_atom_token(lhs)?;
            return Ok(base.set_reaction_type(ReactionFlag::Degradation as u8)
                .set_substrate_key((s as u32) << 16));
        }

        let (s1, s2, s_free) = split_pair(lhs)?;
        let (p1, p2, p_free) = split_pair(rhs)?;
        match (s2, p2) {
            (None, None) => {
                if !same_letter(s1, p1) {
                    return Err(format!("form cannot change in '{rule}'"));
                }
                Ok(base.set_reaction_type(ReactionFlag::Transform as u8)
                    .set_substrate_key((s1 as u32) << 16)
                    .set_product_key((p1 as u32) << 16))
            }
            (Some(s2), Some(p2)) => {
                if !same_letter(s1, p1) || !same_letter(s2, p2) {
                    return Err(format!("form cannot change in '{rule}'"));
                }
                // a bonded pair is one molecule, and a molecule lives in one compartment
                if (!s_free && outside_bit(s1) != outside_bit(s2))
                    || (!p_free && outside_bit(p1) != outside_bit(p2)) {
                    return Err(format!("bonded atoms must be on the same side of the membrane in '{rule}'"));
                }
                let mut t = 0u8;
                if s_free { t |= ReactionFlag::SubstrateContact as u8; }
                if p_free { t |= ReactionFlag::ProductContact as u8; }
                Ok(base.set_reaction_type(t)
                    .set_substrate_key((s1 as u32) << 16 | s2 as u32)
                    .set_product_key((p1 as u32) << 16 | p2 as u32))
            }
            _ => Err(format!("both sides must have the same number of atoms in '{rule}'")),
        }
    }


    pub fn reaction_speed(&self) -> f64 {
        self.reaction_speed
    }

    pub fn id(&self) -> u32 {
        self.id
    }
    pub fn substrate_key(&self) -> u32 {
        self.substrate_key
    }

    pub fn product_key(&self) -> u32 {
        self.product_key
    }
    pub fn set_substrate_key(mut self,substrate_key:u32) -> Self {
        self.substrate_key = substrate_key;
        self
    }
    pub fn set_product_key(mut self, product_key:u32) -> Self {
        self.product_key = product_key;
        self
    }
    pub fn set_reaction_type(mut self, reaction_type: u8) -> Self {
        self.reaction_type = reaction_type;
        self
    }
    pub fn set_reaction_speed(mut self, reaction_speed: f64) -> Self  {
        self.reaction_speed = reaction_speed;
        self
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

    pub fn tags(&self)-> Option<String> {
        self.tags.clone()
    }

    pub fn to_json_value(&self) -> serde_json::Value {
        let mut v = serde_json::json!({
            "s1": format!("0x{:04X}", self.get_substrate_a()),
            "s2": format!("0x{:04X}", self.get_substrate_b()),
            "p1": format!("0x{:04X}", self.get_product_a()),
            "p2": format!("0x{:04X}", self.get_product_b()),
            "reaction_type": format!("0b{:04b}", self.reaction_type),
            "reaction_speed": self.reaction_speed,
        });
        if let Some(t) = &self.tags {
            v["tags"] = serde_json::Value::String(t.clone());
        }
        v
    }
    pub fn to_string(&self) -> String {
        let s1 = key_symbol(self.get_substrate_a());
        let s2 = key_symbol(self.get_substrate_b());
        let p1 = key_symbol(self.get_product_a());
        let p2 = key_symbol(self.get_product_b());
        let sep = |free: bool| if free { "+" } else { "." };

        if self.is_production() {
            format!("0->{s1}")
        } else if self.is_degradation() {
            format!("{s1}->0")
        } else if self.is_transform() {
            format!("{s1}->{p1}")
        } else {
            format!("{s1}{}{s2}->{p1}{}{p2}",
                sep(self.is_substrate_contact()), sep(self.is_product_contact()))
        }
    }
}

/// Parse one atom token like `a1`, `b12` or `A1` into a key (form << 8 | state).
pub fn parse_atom_token(tok: &str) -> Result<u16, String> {
    let mut chars = tok.chars();
    let letter = chars.next().ok_or_else(|| "empty atom".to_string())?;
    // lowercase = inside (a = form 1), uppercase = outside (A = form 1 | OUTSIDE)
    let form = if letter.is_ascii_lowercase() {
        (letter as u8) - 96
    } else if letter.is_ascii_uppercase() {
        ((letter as u8) - 64) | OUTSIDE
    } else {
        return Err(format!("atom '{tok}' must start with a letter"));
    };
    let state: u8 = chars.as_str().parse()
        .map_err(|_| format!("bad state in atom '{tok}'"))?;
    Ok((form as u16) << 8 | state as u16)
}

/// Split one side of a rule: `a1` -> (a1, None, _), `a1+b2` -> (a1, Some(b2), free=true),
/// `a1.b2` -> (a1, Some(b2), free=false).
fn split_pair(side: &str) -> Result<(u16, Option<u16>, bool), String> {
    if let Some((a, b)) = side.split_once('+') {
        Ok((parse_atom_token(a)?, Some(parse_atom_token(b)?), true))
    } else if let Some((a, b)) = side.split_once('.') {
        Ok((parse_atom_token(a)?, Some(parse_atom_token(b)?), false))
    } else {
        Ok((parse_atom_token(side)?, None, false))
    }
}

/// Same letter, ignoring the OUTSIDE (case) bit.
fn same_letter(a: u16, b: u16) -> bool {
    (a >> 8) as u8 & !OUTSIDE == (b >> 8) as u8 & !OUTSIDE
}

fn outside_bit(key: u16) -> bool {
    (key >> 8) as u8 & OUTSIDE != 0
}
