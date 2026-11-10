use std::collections::HashMap;
use std::sync::RwLock;
use once_cell::sync::Lazy;

#[derive(Debug)]
pub struct ReactionMaps {
    forward: HashMap<String, u32>,
    reverse: HashMap<u32, String>,
}

impl ReactionMaps {
    fn new() -> Self {
        let mut forward = HashMap::new();
        let mut reverse = HashMap::new();

        // Your hardcoded mappings
        let reactions = [
            ("terminal", 101),
            ("terminal", 102),
            ("separation", 103),
            ("production", 104),
            ("degradation", 105),
            ("detergent", 106),
            ("detergent unlink", 107),
        ];

        for (name, id) in reactions {
            forward.insert(name.to_string(), id);
            reverse.insert(id, name.to_string());
        }

        Self { forward, reverse }
    }

    pub fn get_id(&self, name: &str) -> Option<u32> {
        self.forward.get(name).copied()
    }

    pub fn get_name(&self, id: u32) -> Option<String> {
        self.reverse.get(&id).cloned()
    }

    pub fn add(&mut self, name: String, id: u32) -> bool {
        if self.forward.contains_key(&name) || self.reverse.contains_key(&id) {
            false
        } else {
            self.forward.insert(name.clone(), id);
            self.reverse.insert(id, name);
            true
        }
    }

    pub fn list_all(&self) -> Vec<(String, u32)> {
        self.forward
            .iter()
            .map(|(k, v)| (k.clone(), *v))
            .collect()
    }
}

// Global singleton
pub static REACTION_MAPS: Lazy<RwLock<ReactionMaps>> =
    Lazy::new(|| RwLock::new(ReactionMaps::new()));

// Convenience functions
pub fn get_reaction_id(name: &str) -> Option<u32> {
    REACTION_MAPS.read().unwrap().get_id(name)
}

pub fn get_reaction_name(id: u32) -> Option<String> {
    REACTION_MAPS.read().unwrap().get_name(id)
}

pub fn add_reaction(name: String, id: u32) -> bool {
    REACTION_MAPS.write().unwrap().add(name, id)
}

pub fn list_reactions() -> Vec<(String, u32)> {
    REACTION_MAPS.read().unwrap().list_all()
}