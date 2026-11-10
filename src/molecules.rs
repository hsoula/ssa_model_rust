use std::collections::{HashMap, HashSet};
use crate::utils::{AtomState,AtomType};
use std::fs::File;
use std::fs;
use std::io::{Write, Result};
use petgraph::graph::Graph;
use petgraph::Undirected;
use petgraph::algo::is_isomorphic_matching;
use serde_json::Value;
use rand::rngs::StdRng;
use rand::Rng;
use crate::trace;

pub type AtomKey = (u8, u8);
pub type MoleculeGraph = Graph<AtomKey, (), Undirected>;
#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: HashMap<u32, (u8, u8)>, // id -> (atom_type, atom_state)
    pub links: HashMap<u32, (u32, u32)>, // link_id -> (atom1_id, atom2_id)
    pub name: String,
    pub formula: String,
}

impl Molecule {
    pub fn new() -> Self {
        Self {
            atoms: HashMap::new(),
            links: HashMap::new(),
            name: String::new(),
            formula: String::new(),
        }
    }

    pub fn add_atom(&mut self, id: u32, atom: (u8, u8)) {
        //let k = self.atoms.len();
        self.atoms.insert(id, atom);
    }

    pub fn size(&self) -> usize {
        self.atoms.len()
    }
    pub fn has_link(&self, link:(u32, u32)) -> bool {
            self.links.iter().any(|(_,l)| l == &link)
    }
    pub fn load_from_file(&mut self, file: &str)  {
        let content = fs::read_to_string(file).unwrap() ;
        self.parse_simple_kcf(content.as_str());
    }
    pub fn get_atoms(&self) -> &HashMap<u32, (u8, u8)> {
        &self.atoms
    }
    pub fn get_links(&self) -> &HashMap<u32, (u32, u32)> {
        &self.links
    }

    pub fn build_graph(&self) -> MoleculeGraph {
        let mut graph = Graph::new_undirected();

        // Créer les nœuds (atomes)
        let mut node_indices = HashMap::new();

        for (index, atom) in self.atoms.iter() {

            let node_idx = graph.add_node(atom.clone() as AtomKey);
            node_indices.insert(index, node_idx);

        }

        // Créer les arêtes (liaisons)
        for (key, bond) in self.links.iter() {
            let node1 = node_indices.get(&bond.0).unwrap();
            let node2 = node_indices.get(&bond.1).unwrap();
                graph.add_edge(*node1, *node2, ());
        }
        graph
    }

    pub fn parse_simple_kcf(&mut self, content: &str)  {
        let lines: Vec<&str> = content.lines().collect();
        let mut i = 0;

        self.atoms.clear();
        self.links.clear();

        while i < lines.len() {
            let line = lines[i].trim();

            if line.is_empty() {
                i += 1;
                continue;
            }

            if line.starts_with("ATOM") {
                self.parse_simple_atoms(&lines, &mut i);
            } else if line.starts_with("BONDS") {
                self.parse_simple_bonds(&lines, &mut i);
            } else if line.starts_with("//") {
                break;
            }

            i += 1;
        }
    }

    fn parse_simple_atoms(&mut self, lines: &[&str], index: &mut usize)  {
        if let Some(count_line) = lines.get(*index) {
            if let Some(count_str) = count_line.strip_prefix("ATOM") {
                let atom_count: usize = count_str.trim().parse().unwrap();

                // Lire les lignes d'atomes
                for i in 1..=atom_count {
                    if let Some(atom_line) = lines.get(*index + i) {
                        let parts: Vec<&str> = atom_line.split_whitespace().collect();
                        if parts.len() >= 3 {
                            let id: u32 = parts[0].parse().unwrap();

                            let type_code: u8 = parts[1].parse().unwrap();
                            let state_code: u8 = parts[2].parse().unwrap();

                            self.atoms.insert(id, (type_code, state_code));
                        } else {
                            panic!("Failed to parse atom line");
                        }
                    } else {
                        panic!("Failed to parse atom line");
                    }
                }
                *index += atom_count;
            }
        }
    }

    fn parse_simple_bonds(&mut self, lines: &[&str], index: &mut usize) {
        // Commencer après la ligne BONDS
        let mut bond_id = 1;
        let mut i = *index + 1;

        while i < lines.len() {
            let line = lines[i].trim();

            if line.is_empty() || line.starts_with("//") {
                break;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                let atom1_id: u32 = parts[0].parse().unwrap();
                let atom2_id: u32 = parts[1].parse().unwrap();

                // Vérifier que les atomes existent
                if !self.atoms.contains_key(&atom1_id) {
                    panic!("Failed to parse atom line");
                }
                if !self.atoms.contains_key(&atom2_id) {
                    panic!("Failed to parse atom line");
                }

                self.links.insert(bond_id, (atom1_id, atom2_id));
                bond_id += 1;
            } else {
                panic!("Failed to parse atom line");
            }

            i += 1;
        }
        *index = i - 1; // Mettre à jour l'index
    }

    pub fn save_to_file(&self, filename: &str) {
        let mut file = File::create(filename).unwrap();

        // Écrire la section ATOM
        writeln!(file, "ATOM {}", self.atoms.len());

        // Trier les atomes par ID pour une sortie
        let mut atom_ids: Vec<&u32> = self.atoms.keys().collect();
        atom_ids.sort();

        for id in atom_ids {
            if let Some(&(val1, val2)) = self.atoms.get(id) {
                writeln!(file, "{} {} {}", id, val1, val2);
            }
        }

        // Écrire la section BONDS
        writeln!(file, "BONDS");

        // Trier les liens par ID
        let mut link_ids: Vec<&u32> = self.links.keys().collect();
        link_ids.sort();

        for (_index, &link_id) in link_ids.iter().enumerate() {
            if let Some(&(atom1, atom2)) = self.links.get(link_id) {
                writeln!(file, "{} {}",atom1, atom2);
            }
        }

        // Marqueur de fin
        writeln!(file, "//");
    }

    // ── JSON loading ─────────────────────────────────────────────────────────────

    /// Build a `Molecule` from a parsed JSON value (KEGG-style format).
    ///
    /// Each atom `{"index": N, "type": "C", "symbol": "C1c", ...}` is mapped:
    ///   - form  = `AtomType::from_kcf_symbol(type)` → u8
    ///   - state = last char of `symbol` parsed as KCF annotation → u8
    pub fn from_json_value(v: &Value) -> Molecule {
        let mut mol = Molecule::new();
        mol.name    = v["entry_id"].as_str().unwrap_or("").to_string();
        mol.formula = v["formula"].as_str().unwrap_or("").to_string();

        if let Some(atoms) = v["atoms"].as_array() {
            for atom in atoms {
                let index  = atom["index"].as_u64().unwrap_or(0) as u32;
                let type_s = atom["type"].as_str().unwrap_or("C");
                let symbol = atom["symbol"].as_str().unwrap_or("");

                let form = AtomType::from_kcf_symbol(type_s).to_u8();

                // Last char of symbol is the KCF state annotation (a/b/c/x/y…)
                let ann = symbol.chars().last()
                    .filter(|c| c.is_ascii_lowercase())
                    .map(|c| {
                        // Return a &str slice pointing at the last byte
                        let n = symbol.len() - c.len_utf8();
                        &symbol[n..]
                    })
                    .unwrap_or("");
                let state = AtomState::from_kcf_annotation(ann).to_u8();

                mol.add_atom(index, (form, state));
            }
        }

        if let Some(bonds) = v["bonds"].as_array() {
            for bond in bonds {
                let index = bond["index"].as_u64().unwrap_or(0) as u32;
                let a1    = bond["atom1"].as_u64().unwrap_or(0) as u32;
                let a2    = bond["atom2"].as_u64().unwrap_or(0) as u32;
                mol.links.insert(index, (a1, a2));
            }
        }

        mol
    }

    /// Load a single molecule from a JSON file.
    pub fn load_from_json_file(path: &str) -> Result<Molecule> {
        let content = fs::read_to_string(path)?;
        let v: Value = serde_json::from_str(&content)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        Ok(Molecule::from_json_value(&v))
    }

    /// Load all `.json` files from a directory, returning one `Molecule` per file.
    /// Files that fail to parse are skipped with a warning on stderr.
    pub fn load_from_json_dir(dir: &str) -> Vec<Molecule> {
        let mut molecules = Vec::new();
        let entries = match fs::read_dir(dir) {
            Ok(e) => e,
            Err(err) => { eprintln!("Cannot read directory {dir}: {err}"); return molecules; }
        };
        for entry in entries.flatten() {
            let p = entry.path();
            if p.extension().and_then(|e| e.to_str()) == Some("json") {
                match Molecule::load_from_json_file(p.to_str().unwrap_or("")) {
                    Ok(mol)  => molecules.push(mol),
                    Err(err) => eprintln!("Skipping {}: {err}", p.display()),
                }
            }
        }
        molecules
    }

    // ── Comparison ───────────────────────────────────────────────────────────────

    /// Count atoms by (form, state) key.
    pub fn atom_composition(&self) -> HashMap<(u8, u8), usize> {
        self.atoms.values()
            .fold(HashMap::new(), |mut acc, &key| { *acc.entry(key).or_insert(0) += 1; acc })
    }

    /// True if this molecule is graph-isomorphic to `other` with matching atom labels.
    /// Two atoms match when they have the same (form, state).
    pub fn is_isomorphic_to(&self, other: &Molecule) -> bool {
        if self.atoms.len() != other.atoms.len() || self.links.len() != other.links.len() {
            return false;
        }
        let g1 = self.build_graph();
        let g2 = other.build_graph();
        is_isomorphic_matching(&g1, &g2, |a, b| a == b, |_, _| true)
    }

    /// Jaccard similarity on atom composition (by type+state count vectors).
    /// Returns 1.0 for identical composition, 0.0 for fully disjoint.
    pub fn composition_similarity(&self, other: &Molecule) -> f64 {
        let c1 = self.atom_composition();
        let c2 = other.atom_composition();

        let all_keys: HashSet<(u8, u8)> = c1.keys().chain(c2.keys()).copied().collect();
        if all_keys.is_empty() { return 1.0; }

        let (intersection, union) = all_keys.iter().fold((0usize, 0usize), |(inter, uni), key| {
            let v1 = *c1.get(key).unwrap_or(&0);
            let v2 = *c2.get(key).unwrap_or(&0);
            (inter + v1.min(v2), uni + v1.max(v2))
        });
        if union == 0 { 1.0 } else { intersection as f64 / union as f64 }
    }

    /// Compare a list of molecules against a reference and group by isomorphism class.
    /// Returns a list of (representative, members) pairs — one per distinct class.
    pub fn group_by_isomorphism(molecules: &[Molecule]) -> Vec<(&Molecule, Vec<usize>)> {
        let mut classes: Vec<(&Molecule, Vec<usize>)> = Vec::new();
        'outer: for (i, mol) in molecules.iter().enumerate() {
            for (rep, members) in &mut classes {
                if mol.is_isomorphic_to(rep) {
                    members.push(i);
                    continue 'outer;
                }
            }
            classes.push((mol, vec![i]));
        }
        classes
    }

    // ── Distance metrics ─────────────────────────────────────────────────────────

    /// 1 − composition_similarity (Jaccard on atom type+state counts).
    pub fn composition_distance(&self, other: &Molecule) -> f64 {
        1.0 - self.composition_similarity(other)
    }

    /// Normalised L1 distance on sorted node-degree sequences.
    /// Both sequences are zero-padded to equal length before comparison.
    pub fn degree_sequence_distance(&self, other: &Molecule) -> f64 {
        let deg = |mol: &Molecule| -> Vec<usize> {
            let mut d: HashMap<u32, usize> = mol.atoms.keys().map(|&k| (k, 0)).collect();
            for &(a, b) in mol.links.values() {
                *d.entry(a).or_insert(0) += 1;
                *d.entry(b).or_insert(0) += 1;
            }
            let mut v: Vec<usize> = d.values().copied().collect();
            v.sort_unstable();
            v
        };

        let mut d1 = deg(self);
        let mut d2 = deg(other);

        // zero-pad shorter sequence
        let max_len = d1.len().max(d2.len());
        d1.resize(max_len, 0);
        d2.resize(max_len, 0);

        let l1: usize = d1.iter().zip(d2.iter()).map(|(&a, &b)| a.abs_diff(b)).sum();
        let norm: usize = d1.iter().chain(d2.iter()).copied().sum::<usize>().max(1);
        l1 as f64 / norm as f64
    }

    /// Weisfeiler-Lehman graph distance (1 − Jaccard of WL colour histograms).
    /// `depth` controls the number of refinement rounds (1–3 is usually enough).
    pub fn wl_distance(&self, other: &Molecule, depth: usize) -> f64 {
        let h1 = wl_histogram(self, depth);
        let h2 = wl_histogram(other, depth);

        let all_keys: HashSet<&String> = h1.keys().chain(h2.keys()).collect();
        if all_keys.is_empty() { return 0.0; }

        let (inter, union) = all_keys.iter().fold((0usize, 0usize), |(i, u), k| {
            let v1 = *h1.get(*k).unwrap_or(&0);
            let v2 = *h2.get(*k).unwrap_or(&0);
            (i + v1.min(v2), u + v1.max(v2))
        });
        if union == 0 { 0.0 } else { 1.0 - inter as f64 / union as f64 }
    }

    // ── Random molecule ───────────────────────────────────────────────────────────

    /// Generate a random connected molecule with `n_atoms` atoms.
    ///
    /// - Atom forms are sampled uniformly from `forms`; states are always 0.
    /// - Connectivity is a random spanning tree (Prüfer sequence) → always connected.
    /// - `extra_edges` additional random edges are added (may create cycles).
    pub fn random(n_atoms: usize, forms: &[u8], states: &[u8], extra_edges: usize, rng: &mut StdRng) -> Molecule {
        let mut mol = Molecule::new();
        if n_atoms == 0 || forms.is_empty() { return mol; }

        for id in 1..=(n_atoms as u32) {
            let form = forms[rng.random_range(0..forms.len())];
            let state = states[rng.random_range(0..states.len())];
            mol.add_atom(id, (form, state));
        }

        // Random spanning tree: connect each new node to a random existing node.
        // This always produces a connected tree with exactly n_atoms−1 edges.
        if n_atoms > 1 {
            let mut bond_id = 1u32;
            for new_node in 2..=(n_atoms as u32) {
                let existing = rng.random_range(1..new_node);
                mol.links.insert(bond_id, (existing, new_node));
                bond_id += 1;
            }

            // Extra random edges (no self-loops, no duplicates)
            let tree_edges: HashSet<(u32, u32)> = mol.links.values()
                .map(|&(a, b)| (a.min(b), a.max(b)))
                .collect();
            let mut used = tree_edges;
            let mut tried = 0usize;
            let mut added = 0usize;
            while added < extra_edges && tried < extra_edges * 10 {
                let a = rng.random_range(1..=(n_atoms as u32));
                let b = rng.random_range(1..=(n_atoms as u32));
                tried += 1;
                if a == b { continue; }
                let key = (a.min(b), a.max(b));
                if used.contains(&key) { continue; }
                used.insert(key);
                mol.links.insert(bond_id, (a, b));
                bond_id += 1;
                added += 1;
            }
        }

        mol.name    = format!("rand_{}", n_atoms);
        mol.formula = String::new();
        mol
    }
}

// ── WL histogram (free function) ─────────────────────────────────────────────

/// Compute a Weisfeiler-Lehman colour histogram for `mol` up to `depth` rounds.
fn wl_histogram(mol: &Molecule, depth: usize) -> HashMap<String, usize> {
    let n = mol.atoms.len();
    if n == 0 { return HashMap::new(); }

    // Adjacency: atom_id → neighbour ids
    let mut adj: HashMap<u32, Vec<u32>> = mol.atoms.keys().map(|&k| (k, Vec::new())).collect();
    for &(a, b) in mol.links.values() {
        adj.entry(a).or_default().push(b);
        adj.entry(b).or_default().push(a);
    }

    // Initial colour = "form_state"
    let mut colours: HashMap<u32, String> = mol.atoms.iter()
        .map(|(&id, &(f, s))| (id, format!("{}_{}", f, s)))
        .collect();

    let mut hist: HashMap<String, usize> = HashMap::new();
    for c in colours.values() { *hist.entry(c.clone()).or_insert(0) += 1; }

    for _ in 0..depth {
        let mut next: HashMap<u32, String> = HashMap::new();
        for (&id, nbrs) in &adj {
            let mut nbr_colours: Vec<&str> = nbrs.iter()
                .map(|n| colours.get(n).map(|s| s.as_str()).unwrap_or(""))
                .collect();
            nbr_colours.sort_unstable();
            let new_col = format!("{}|{}", colours[&id], nbr_colours.join(","));
            *hist.entry(new_col.clone()).or_insert(0) += 1;
            next.insert(id, new_col);
        }
        colours = next;
    }
    hist
}