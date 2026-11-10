use std::iter::Iterator;
use std::string::{String, ToString};
use std::vec::Vec;
use rand::Rng;
use rand::rngs::StdRng;
use crate::Reaction;

/*
  Gene8 alphabet — 16 symbols (one hex digit each, 0x00–0x0F)
  ──────────────────────────────────────────────────────────────
  E (0x0E) = START marker   — opens a reaction pattern
  F (0x0F) = END   marker   — closes a reaction pattern

  Payload tokens (appear between E…F):
    0,1  →  sf1++   substrate 1 form count  (×2 probability)
    2,3  →  ss1++   substrate 1 state count (×2 probability)
    4,5  →  sf2++   substrate 2 form count  (×2 probability)
    6,7  →  ss2++   substrate 2 state count (×2 probability)
    8    →  contact++  reaction-type bits (mod 4)
    9,A  →  ps1++   product 1 state count  (×2 probability, form inherited from sf1)
    B,C  →  ps2++   product 2 state count  (×2 probability, form inherited from sf2)
    D    →  speed += 0.1

  Product forms are always inherited from substrate forms:
    p1.form = sf1,  p2.form = sf2
  Only states can differ between substrate and product.

  Form counters (sf1, sf2) are capped at `max_count` (default 10 — matches the
  10 named AtomTypes: C=1 … I=10).  State counters (ss1, ss2, ps1, ps2) are
  capped at `max_count × 2`, giving twice the range (states are assumed more
  numerous than atom types).

  A reaction is only valid if both sf1 > 0 and sf2 > 0.
 */

#[derive(Clone)]
pub struct Gene8 {
    pub genome: Vec<u8>,
    pub max_count: u8,
}

impl Gene8 {
    pub fn new(size: usize, rng: &mut StdRng) -> Self {
        Self::new_with_max(size, 10, rng)
    }

    pub fn new_with_max(size: usize, max_count: u8, rng: &mut StdRng) -> Self {
        let mut genome = vec![0u8; size];
        rng.fill(genome.as_mut_slice());
        genome = genome.iter().map(|x| x % 16).collect();
        Self { genome, max_count }
    }

    pub fn len(&self) -> usize { self.genome.len() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn to_string(&self) -> String {
        self.genome.iter().map(|x| format!("{:X}", x)).collect()
    }

    pub fn extract_patterns_efficient(&self) -> Vec<&[u8]> {
        let mut patterns = Vec::new();
        for (i, &byte) in self.genome.iter().enumerate() {
            if byte == 0x0E {
                'out: for j in (i + 1)..self.genome.len() {
                    if self.genome[j] == 0x0F {
                        patterns.push(&self.genome[i..=j]);
                        break 'out;
                    }
                }
            }
        }
        patterns
    }

    pub fn point_mutations(&self, point_mutation_prob: f64, rng: &mut StdRng) -> Gene8 {
        let fp = self.genome.iter().map(|x|
            if rng.random_bool(point_mutation_prob) {
                rng.random_range(0u8..16u8)
            } else {
                *x
            }).collect();
        Gene8 { genome: fp, max_count: self.max_count }
    }

    pub fn point_addition(&self, point_mutation_prob: f64, rng: &mut StdRng) -> Gene8 {
        let mut gene_muted = self.genome.clone();
        let n = gene_muted.len();
        let r = ((n as f64) * rng.random_range(0f64..1f64)) as u32;
        for _i in 0..r {
            let u = gene_muted.len();
            let k = rng.random_range(0usize..u);
            let v = rng.random_range(0u8..16u8);
            gene_muted.insert(k, v);
        }
        Gene8 { genome: gene_muted, max_count: self.max_count }
    }

    pub fn to_reactions(&self) -> Vec<crate::reactions::Reaction> {
        self.extract_patterns_efficient()
            .iter()
            .enumerate()
            .map(|(i, pattern)| {
                crate::reactions::Reaction::new_from_vec8(i as u32, pattern.to_vec(), self.max_count)
            })
            .filter(|r| {
                (r.get_substrate_a() >> 8) > 0 && (r.get_substrate_b() >> 8) > 0 && r.substrate_key() != r.product_key()
            })
            .collect()
    }

    pub fn save_to_json(&self, path: &str) {
        let reactions = self.to_reactions();
        let json: Vec<serde_json::Value> = reactions.iter()
            .map(|r| r.to_json_value())
            .collect();
        std::fs::write(path, serde_json::to_string_pretty(&json).unwrap()).unwrap();
    }

    pub fn crossing(&self, other: &Gene8, rng: &mut StdRng) -> (Gene8, Gene8) {
        let n1 = self.len();
        let n2 = other.len();
        let u1 = ((n1 as f64) * rng.random_range(0f64..1f64)) as usize;
        let u2 = ((n2 as f64) * rng.random_range(0f64..1f64)) as usize;
        let mut v1 = self.genome[0..u1].to_vec();
        let mut v2 = other.genome[0..u2].to_vec();
        v1.extend_from_slice(&other.genome[u2..n2]);
        v2.extend_from_slice(&self.genome[u1..n1]);
        (Gene8 { genome: v1, max_count: self.max_count },
         Gene8 { genome: v2, max_count: self.max_count })
    }
}
