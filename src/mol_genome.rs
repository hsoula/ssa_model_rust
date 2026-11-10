//! MolGenome — genome encoding for molecules.
//!
//! # Alphabet (16 nibble tokens, 0x00–0x0F)
//!
//! | Token | Symbol | Role |
//! |-------|--------|------|
//! | 0x0E  | `E`    | START — opens a molecule pattern |
//! | 0x0F  | `F`    | END   — closes a molecule pattern |
//! | 0x00  | `0`    | State token — increments state of the most recent atom |
//! | 0x01–0x0D | `1`–`D` | Form tokens — creates a new atom of that form |
//!
//! # State encoding
//!
//! The state of each atom equals **(number of 0-tokens that immediately follow
//! the form token, until the next form token or `F`) − 1**, clamped to
//! `[0, max_state]`.
//!
//! In other words: at least one `0` is needed to represent state 0; each
//! additional `0` adds 1 to the state.  If no `0` follows a form token,
//! the atom gets state 0.
//!
//! Example: `E100201000F`
//!   inner = `1 0 0 2 0 1 0 0 0`
//!   → atom1: form=1('a'), 2 trailing zeros → state=1 → `a1`
//!   → atom2: form=2('b'), 1 trailing zero  → state=0 → `b0`
//!   → atom3: form=1('a'), 3 trailing zeros → state=2 → `a2`
//!   bonds: a1—b0—a2 (linear chain)
//!
//! # Bond structure
//!
//! Atoms within a pattern are bonded as a **linear chain**:
//! atom[0]—atom[1]—atom[2]—…
//!
//! # Form clamping
//!
//! A form token value `t ∈ [1,13]` is mapped to a form in `[1, max_form]` by
//! `((t − 1) % max_form) + 1`.  All non-zero tokens therefore contribute an atom.

use rand::Rng;
use rand::rngs::StdRng;
use crate::molecules::Molecule;

const START: u8 = 0x0E;
const END:   u8 = 0x0F;
const STATE: u8 = 0x00;

#[derive(Clone, Debug)]
pub struct MolGenome {
    /// Nibble values 0x00–0x0F.
    pub genome: Vec<u8>,
    /// Maximum form index (forms are 1..=max_form; token values are wrapped).
    pub max_form: u8,
    /// Maximum state value (state is clamped to [0, max_state]).
    pub max_state: u8,
}

impl MolGenome {
    // ── Constructors ──────────────────────────────────────────────────────────

    /// Random genome of `size` nibble tokens.
    pub fn new(size: usize, max_form: u8, max_state: u8, rng: &mut StdRng) -> Self {
        let genome: Vec<u8> = (0..size).map(|_| rng.random_range(0u8..16u8)).collect();
        Self { genome, max_form, max_state }
    }

    /// Empty genome (useful for hand-crafting test cases).
    pub fn empty(max_form: u8, max_state: u8) -> Self {
        Self { genome: Vec::new(), max_form, max_state }
    }

    /// Build a genome from a hex string (e.g. `"E100201000F"`).
    ///
    /// Unknown characters are silently skipped.
    pub fn from_hex(s: &str, max_form: u8, max_state: u8) -> Self {
        let genome: Vec<u8> = s.chars()
            .filter_map(|c| u8::from_str_radix(&c.to_string(), 16).ok())
            .collect();
        Self { genome, max_form, max_state }
    }

    // ── Display ───────────────────────────────────────────────────────────────

    /// Hexadecimal string representation (e.g. `"E100201000F"`).
    pub fn to_hex(&self) -> String {
        self.genome.iter().map(|x| format!("{:X}", x)).collect()
    }

    /// Annotated view: each token on its own line with its role.
    pub fn annotate(&self) -> String {
        let mut out = String::new();
        for &b in &self.genome {
            let label = match b {
                START => "[START]",
                END   => "[END]",
                STATE => "[state+1]",
                _     => "[form]",
            };
            out.push_str(&format!("{:X} {label}\n", b));
        }
        out
    }

    // ── Pattern extraction ────────────────────────────────────────────────────

    /// Return all `E…F` sub-slices (inclusive of the E and F markers).
    pub fn extract_patterns(&self) -> Vec<&[u8]> {
        let mut patterns = Vec::new();
        let g = &self.genome;
        for (i, &b) in g.iter().enumerate() {
            if b == START {
                for j in (i + 1)..g.len() {
                    if g[j] == END {
                        patterns.push(&g[i..=j]);
                        break;
                    }
                }
            }
        }
        patterns
    }

    // ── Decoding ──────────────────────────────────────────────────────────────

    /// Decode a single `E…F` pattern (including delimiters) into a `Molecule`.
    ///
    /// Returns `None` if the pattern contains fewer than 2 atoms.
    pub fn decode_pattern(pattern: &[u8], max_form: u8, max_state: u8) -> Option<Molecule> {
        // Strip delimiters
        if pattern.len() < 2 { return None; }
        let inner: &[u8] = if pattern[0] == START && pattern[pattern.len()-1] == END {
            &pattern[1..pattern.len()-1]
        } else {
            pattern
        };

        // Locate form token positions
        let form_positions: Vec<usize> = inner.iter().enumerate()
            .filter(|&(_, b)| *b != STATE && *b != START && *b != END)
            .map(|(i, _)| i)
            .collect();

        if form_positions.len() < 2 { return None; }

        let mut mol = Molecule::new();

        for (atom_idx, &fp) in form_positions.iter().enumerate() {
            let raw_form = inner[fp];
            // Clamp: map form tokens 1..=0x0D into [1, max_form] via wrap
            let form = ((raw_form - 1) % max_form) + 1;

            // Count 0-tokens immediately following until next form token or end
            let next_form_pos = form_positions.get(atom_idx + 1).copied()
                .unwrap_or(inner.len());
            let n_zeros = inner[fp+1..next_form_pos]
                .iter()
                .filter(|&&b| b == STATE)
                .count();

            // state = n_zeros - 1 (min 0, capped at max_state)
            let state = if n_zeros == 0 { 0u8 }
                        else { ((n_zeros as u8) - 1).min(max_state) };

            let atom_id = atom_idx as u32 + 1;
            mol.atoms.insert(atom_id, (form, state));

            // Linear bond to previous atom
            if atom_idx > 0 {
                mol.links.insert(atom_idx as u32 - 1, (atom_id - 1, atom_id));
            }
        }

        mol.name = String::from("mol");
        Some(mol)
    }

    /// Decode all valid patterns to molecules (requires ≥ 2 atoms each).
    pub fn to_molecules(&self) -> Vec<Molecule> {
        self.extract_patterns()
            .iter()
            .enumerate()
            .filter_map(|(i, &pat)| {
                let mut m = Self::decode_pattern(pat, self.max_form, self.max_state)?;
                m.name = format!("mol_{i}");
                Some(m)
            })
            .collect()
    }

    // ── Mutation & crossover ──────────────────────────────────────────────────

    /// Point mutation: each token independently replaced with prob `rate`.
    pub fn mutate(&self, rate: f64, rng: &mut StdRng) -> MolGenome {
        let genome = self.genome.iter().map(|&b| {
            if rng.random_bool(rate) { rng.random_range(0u8..16u8) } else { b }
        }).collect();
        MolGenome { genome, max_form: self.max_form, max_state: self.max_state }
    }

    /// Single-point crossover: returns two children.
    pub fn cross(&self, other: &MolGenome, rng: &mut StdRng) -> (MolGenome, MolGenome) {
        let n1 = self.genome.len();
        let n2 = other.genome.len();
        let cut1 = rng.random_range(0..=n1);
        let cut2 = rng.random_range(0..=n2);

        let mut g1 = self.genome[..cut1].to_vec();
        g1.extend_from_slice(&other.genome[cut2..]);

        let mut g2 = other.genome[..cut2].to_vec();
        g2.extend_from_slice(&self.genome[cut1..]);

        (MolGenome { genome: g1, max_form: self.max_form, max_state: self.max_state },
         MolGenome { genome: g2, max_form: self.max_form, max_state: self.max_state })
    }

    /// Insertion mutation: insert `rate * len` random tokens at random positions.
    pub fn insert(&self, rate: f64, rng: &mut StdRng) -> MolGenome {
        let mut g = self.genome.clone();
        let n_insert = ((g.len() as f64 * rate * rng.random::<f64>()) as usize).max(0);
        for _ in 0..n_insert {
            let pos = rng.random_range(0..=g.len());
            let tok = rng.random_range(0u8..16u8);
            g.insert(pos, tok);
        }
        MolGenome { genome: g, max_form: self.max_form, max_state: self.max_state }
    }

    /// Deletion mutation: delete `rate * len` random tokens.
    pub fn delete(&self, rate: f64, rng: &mut StdRng) -> MolGenome {
        if self.genome.is_empty() { return self.clone(); }
        let mut g = self.genome.clone();
        let n_del = ((g.len() as f64 * rate * rng.random::<f64>()) as usize).max(0);
        for _ in 0..n_del {
            if g.is_empty() { break; }
            let pos = rng.random_range(0..g.len());
            g.remove(pos);
        }
        MolGenome { genome: g, max_form: self.max_form, max_state: self.max_state }
    }

    // ── Utility ───────────────────────────────────────────────────────────────

    pub fn len(&self) -> usize { self.genome.len() }
    pub fn is_empty(&self) -> bool { self.genome.is_empty() }

    /// Count valid molecule patterns in this genome.
    pub fn n_patterns(&self) -> usize { self.extract_patterns().len() }

    /// Count molecules that decode to ≥ 2 atoms.
    pub fn n_valid_molecules(&self) -> usize { self.to_molecules().len() }

    /// Serialize to a JSON value: `{ genome, max_form, max_state, molecules: [...] }`.
    pub fn to_json(&self) -> serde_json::Value {
        let mols: Vec<serde_json::Value> = self.to_molecules().iter()
            .map(|m| {
                let mut atoms: Vec<serde_json::Value> = m.atoms.iter()
                    .map(|(&id, &(f, s))| serde_json::json!({"id": id, "form": f, "state": s}))
                    .collect();
                atoms.sort_by_key(|v| v["id"].as_u64().unwrap_or(0));
                let mut bonds: Vec<serde_json::Value> = m.links.iter()
                    .map(|(&id, &(a, b))| serde_json::json!({"id": id, "a": a, "b": b}))
                    .collect();
                bonds.sort_by_key(|v| v["id"].as_u64().unwrap_or(0));
                serde_json::json!({
                    "name": m.name, "n_atoms": m.atoms.len(), "n_bonds": m.links.len(),
                    "atoms": atoms, "bonds": bonds,
                })
            })
            .collect();
        serde_json::json!({
            "genome":     self.to_hex(),
            "genome_len": self.genome.len(),
            "max_form":   self.max_form,
            "max_state":  self.max_state,
            "n_patterns": self.n_patterns(),
            "n_molecules": mols.len(),
            "molecules":  mols,
        })
    }

    /// Save genome + decoded molecules to a JSON file.
    pub fn save_to_json(&self, path: &str) {
        std::fs::write(path, serde_json::to_string_pretty(&self.to_json()).unwrap())
            .unwrap_or_else(|e| eprintln!("Cannot write {path}: {e}"));
    }

    /// Render a human-readable summary of all decoded molecules.
    pub fn describe(&self) -> String {
        let mols = self.to_molecules();
        if mols.is_empty() {
            return format!("genome {} | no valid molecules", self.to_hex());
        }
        let descs: Vec<String> = mols.iter().map(|m| {
            let mut atoms: Vec<(u32, u8, u8)> = m.atoms.iter()
                .map(|(&id, &(f, s))| (id, f, s)).collect();
            atoms.sort_by_key(|(id, _, _)| *id);
            let atom_str = atoms.iter()
                .map(|(_, f, s)| format!("{}{}", (f + 96) as char, s))
                .collect::<Vec<_>>().join("-");
            format!("[{}]", atom_str)
        }).collect();
        format!("genome {} | {} mol(s): {}",
            self.to_hex(), mols.len(), descs.join(", "))
    }
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify the canonical example from the module doc:
    /// E100201000F → a1-b0-a2
    #[test]
    fn test_canonical_example() {
        let g = MolGenome::from_hex("E100201000F", 2, 4);
        let mols = g.to_molecules();
        assert_eq!(mols.len(), 1, "expected 1 molecule");
        let m = &mols[0];
        assert_eq!(m.atoms.len(), 3);
        assert_eq!(m.links.len(), 2);
        // atom 1: form=1, state=1  (a1)
        assert_eq!(m.atoms[&1], (1, 1));
        // atom 2: form=2, state=0  (b0)
        assert_eq!(m.atoms[&2], (2, 0));
        // atom 3: form=1, state=2  (a2)
        assert_eq!(m.atoms[&3], (1, 2));
        // bonds: 1-2 and 2-3
        let bond_set: std::collections::HashSet<(u32,u32)> = m.links.values().cloned().collect();
        assert!(bond_set.contains(&(1,2)));
        assert!(bond_set.contains(&(2,3)));
    }

    #[test]
    fn test_single_atom_pattern_ignored() {
        // Only one atom between E…F → not a valid molecule
        let g = MolGenome::from_hex("E10F", 2, 4);
        assert_eq!(g.to_molecules().len(), 0);
    }

    #[test]
    fn test_no_zeros_state_is_zero() {
        // E12F → a0-b0 (no zeros, so state defaults to 0)
        let g = MolGenome::from_hex("E12F", 2, 4);
        let mols = g.to_molecules();
        assert_eq!(mols.len(), 1);
        let m = &mols[0];
        assert_eq!(m.atoms[&1], (1, 0)); // a0
        assert_eq!(m.atoms[&2], (2, 0)); // b0
    }

    #[test]
    fn test_state_clamped_to_max() {
        // E100000F: form=1, then 5 zeros → state = 5-1 = 4
        // But if max_state=2, it clamps to 2
        let g = MolGenome::from_hex("E1000002F", 2, 2);
        let mols = g.to_molecules();
        assert_eq!(mols.len(), 1);
        let m = &mols[0];
        assert_eq!(m.atoms[&1].1, 2); // clamped at max_state=2
    }

    #[test]
    fn test_form_wrapping() {
        // max_form=2; token 3 → ((3-1) % 2) + 1 = 1 = form 'a'
        // E30F: tokens 3 and 0 → atom form=1, state=0; then nothing else
        // Not valid (only 1 atom), but test wrapping via 2-atom genome
        let g = MolGenome::from_hex("E310F", 2, 4);
        let mols = g.to_molecules();
        assert_eq!(mols.len(), 1);
        // token 3 wraps to form ((3-1)%2)+1 = 1
        assert_eq!(mols[0].atoms.get(&1).map(|&(f, _)| f), Some(1));
    }

    #[test]
    fn test_multiple_patterns() {
        // Two patterns: E10F and E200F
        let g = MolGenome::from_hex("E12FE200F", 2, 4);
        // First pattern E12F: a0-b0 (valid)
        // Second pattern E200F: form=2(b), then nothing for state (0 zeros), form=2(b), one zero for state → b0-b0
        // Wait: E200F inner = 2 0 0 → two form positions: token '2' at pos 0, token '0' is STATE, token '0' is STATE
        // But '2' is the only non-state token... actually inner='200' in hex:
        //   0x02 = form token, 0x00 = state, 0x00 = state
        // Only 1 form token → not a valid molecule
        // Hmm, let me think about E2020F:
        //   inner = 0x02 0x00 0x02 0x00 → form=2 then 1 zero then form=2 then 1 zero
        //   atom1: form=2, 1 zero → state=0 → b0
        //   atom2: form=2, 1 zero → state=0 → b0
        // So E12FE2020F:
        let g2 = MolGenome::from_hex("E12FE2020F", 2, 4);
        let mols = g2.to_molecules();
        assert_eq!(mols.len(), 2);
    }

    #[test]
    fn test_describe() {
        let g = MolGenome::from_hex("E100201000F", 2, 4);
        let desc = g.describe();
        assert!(desc.contains("a1-b0-a2"), "describe should show a1-b0-a2, got: {desc}");
    }

    #[test]
    fn test_mutate_preserves_length() {
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(0);
        let g = MolGenome::new(50, 3, 4, &mut rng);
        let m = g.mutate(0.1, &mut rng);
        assert_eq!(m.len(), g.len());
    }

    #[test]
    fn test_cross_lengths() {
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(1);
        let g1 = MolGenome::new(30, 2, 2, &mut rng);
        let g2 = MolGenome::new(40, 2, 2, &mut rng);
        let (c1, c2) = g1.cross(&g2, &mut rng);
        assert_eq!(c1.len() + c2.len(), g1.len() + g2.len(),
            "crossover should preserve total token count");
    }
}
