use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use serde_with::{serde_as, DisplayFromStr};

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
    pub reaction_speed: f64
}
