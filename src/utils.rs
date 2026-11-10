use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use serde_with::{serde_as, DisplayFromStr};
use crate::trace;
use std::sync::atomic::{AtomicU64, Ordering};


fn serialize_hex<S>(value: &u16, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_str(&format!("0x{:04X}", value))
}

fn deserialize_hex<'de, D>(deserializer: D) -> Result<u16, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;

    if let Some(hex_str) = s.strip_prefix("0x") {
        u16::from_str_radix(hex_str, 16)
            .map_err(serde::de::Error::custom)
    } else {
        s.parse::<u16>()
            .map_err(serde::de::Error::custom)
    }
}

fn serialize_bin<S>(value: &u8, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_str(&format!("0b{:04b}", value))
}

fn deserialize_bin<'de, D>(deserializer: D) -> Result<u8, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;

    if let Some(bin_str) = s.strip_prefix("0b") {
        u8::from_str_radix(bin_str, 2)
            .map_err(serde::de::Error::custom)
    } else {
        s.parse::<u8>()
            .map_err(serde::de::Error::custom)
    }
}

#[serde_as]
#[derive(Debug, Deserialize)]
pub struct JsonReaction {
    #[serde(serialize_with = "serialize_hex", deserialize_with = "deserialize_hex")]
    pub s1: u16,
    #[serde(serialize_with = "serialize_hex", deserialize_with = "deserialize_hex")]
    pub s2: u16,
    #[serde(serialize_with = "serialize_hex", deserialize_with = "deserialize_hex")]
    pub p1: u16,
    #[serde(serialize_with = "serialize_hex", deserialize_with = "deserialize_hex")]
    pub p2: u16,
    #[serde(serialize_with = "serialize_bin", deserialize_with = "deserialize_bin")]
    pub reaction_type: u8,
    pub reaction_speed: f64,

    pub tags: Option<String>,
}



/// Parse a `KEY:COUNT` atom spec into one or more `(key, count)` pairs.
///
/// `KEY` is normally a 4-hex-digit `0xFFSS` byte pair (form, state), matching
/// the packed atom key (`form << 8 | state`). Either byte's 2-digit group may
/// contain a `*` to wildcard that byte: a wildcarded form expands to
/// `1..=max_form`, a wildcarded state expands to `0..=max_state`, and the
/// result is the cross product (so `0x0*0*` yields `max_form * max_state`
/// pairs, all sharing `COUNT`).
///
/// Examples: `0x0100:300` (form=1,state=0 → 300 atoms), `0x0*03:50` (every
/// form, state=3, 50 atoms each), `0x01**:50` (form=1, every state).
/// A plain decimal `KEY:COUNT` (no `0x` prefix) is also accepted, unchanged.
pub fn parse_atom_spec(s: &str, max_form: u8, max_state: u8) -> Result<Vec<(u16, usize)>, String> {
    let (key_str, count_str) = s.split_once(':')
        .ok_or_else(|| format!("expected KEY:COUNT, got '{s}'"))?;
    let count = count_str.trim().parse::<usize>().map_err(|e| e.to_string())?;
    let key_str = key_str.trim();

    let hex = match key_str.strip_prefix("0x").or_else(|| key_str.strip_prefix("0X")) {
        Some(hex) => hex,
        None => {
            let key = key_str.parse::<u16>().map_err(|e| e.to_string())?;
            return Ok(vec![(key, count)]);
        }
    };

    if !hex.contains('*') {
        let key = u16::from_str_radix(hex, 16).map_err(|e| e.to_string())?;
        return Ok(vec![(key, count)]);
    }

    if hex.len() != 4 {
        return Err(format!(
            "wildcard key '{key_str}' must have exactly 4 hex digits (e.g. 0x0*0*), got '{hex}'"
        ));
    }
    let (form_part, state_part) = hex.split_at(2);

    let forms: Vec<u8> = if form_part.contains('*') {
        (1..=max_form).collect()
    } else {
        vec![u8::from_str_radix(form_part, 16).map_err(|e| e.to_string())?]
    };
    let states: Vec<u8> = if state_part.contains('*') {
        (0..=max_state).collect()
    } else {
        vec![u8::from_str_radix(state_part, 16).map_err(|e| e.to_string())?]
    };

    Ok(forms.iter()
        .flat_map(|&f| states.iter().map(move |&st| (((f as u16) << 8) | st as u16, count)))
        .collect())
}

// Générateur d'IDs global et thread-safe
static NEXT_ID: AtomicU64 = AtomicU64::new(1);

pub type ObjectId = u64;

pub fn generate_id() -> ObjectId {
    NEXT_ID.fetch_add(1, Ordering::Relaxed)
}


#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomType {
    C,  // Carbon
    N,  // Azote
    O,  // Oxygen
    P,  // Phosphore
    S,  // Soufre
    H,  // Hydrogen
    F,  // Fluor
    Cl, // Chlore
    Br, // Brome
    I,  // Iode
    Other(u8), // other
}

impl AtomType {
    pub fn from_kcf_symbol(symbol: &str) -> Self {
        match symbol {
            "C" => AtomType::C,
            "N" => AtomType::N,
            "O" => AtomType::O,
            "P" => AtomType::P,
            "S" => AtomType::S,
            "H" => AtomType::H,
            "F" => AtomType::F,
            "Cl" => AtomType::Cl,
            "Br" => AtomType::Br,
            "I" => AtomType::I,
            _ => {
                // Pour les symboles complexes comme "C8y", "N4y", etc.
                // On prend le premier caractère comme type d'atome
                if let Some(first_char) = symbol.chars().next() {
                    match first_char {
                        'C' => AtomType::C,
                        'N' => AtomType::N,
                        'O' => AtomType::O,
                        'P' => AtomType::P,
                        'S' => AtomType::S,
                        'H' => AtomType::H,
                        'F' => AtomType::F,
                        'C' if symbol.starts_with("Cl") => AtomType::Cl,
                        'B' if symbol.starts_with("Br") => AtomType::Br,
                        'I' => AtomType::I,
                        _ => AtomType::Other(0),
                    }
                } else {
                    AtomType::Other(0)
                }
            }
        }
    }

    pub fn to_u8(self) -> u8 {
        match self {
            AtomType::C => 1,
            AtomType::N => 2,
            AtomType::O => 3,
            AtomType::P => 4,
            AtomType::S => 5,
            AtomType::H => 6,
            AtomType::F => 7,
            AtomType::Cl => 8,
            AtomType::Br => 9,
            AtomType::I => 10,
            AtomType::Other(code) => code,
        }
    }

    pub fn from_u8(code: u8) -> Self {
        match code {
            1 => AtomType::C,
            2 => AtomType::N,
            3 => AtomType::O,
            4 => AtomType::P,
            5 => AtomType::S,
            6 => AtomType::H,
            7 => AtomType::F,
            8 => AtomType::Cl,
            9 => AtomType::Br,
            10 => AtomType::I,
            n => AtomType::Other(n),
        }
    }
}

// État de l'atome selon KCF
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AtomState {
    Normal,      // État normal
    Aromatic,    // Aromatique
    Charged,     // Chargé
    Radical,     // Radical
    Isotope,     // Isotope
    Other(u8),   // Autres états
}

impl AtomState {
    pub fn from_kcf_annotation(annotation: &str) -> Self {
        match annotation {
            "" => AtomState::Normal,
            "y" => AtomState::Aromatic,  // Exemple: C8y
            "x" => AtomState::Charged,   // Exemple: N5x
            "a" => AtomState::Radical,   // Exemple: O1a
            "b" => AtomState::Isotope,   // Exemple: P1b
            "c" => AtomState::Other(1),  // Exemple: O1c
            _ => AtomState::Other(2),
        }
    }

    pub fn to_u8(self) -> u8 {
        match self {
            AtomState::Normal => 0,
            AtomState::Aromatic => 1,
            AtomState::Charged => 2,
            AtomState::Radical => 3,
            AtomState::Isotope => 4,
            AtomState::Other(code) => code,
        }
    }

    pub fn from_u8(code: u8) -> Self {
        match code {
            0 => AtomState::Normal,
            1 => AtomState::Aromatic,
            2 => AtomState::Charged,
            3 => AtomState::Radical,
            4 => AtomState::Isotope,
            n => AtomState::Other(n),
        }
    }
}
