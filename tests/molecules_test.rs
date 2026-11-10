use ssa_model::utils::{AtomType, AtomState};
use ssa_model::molecules::{Molecule};

#[cfg(test)]
mod tests {
    use ssa_model::trace;
    use super::*;

    #[test]
    fn test_parse_file() {
        let mut molecule = Molecule::new();
        molecule.load_from_file("data/reactor1.mol");
        let x = molecule.get_atoms().iter().count();
        assert_eq!(x, 5);
        let x = molecule.get_links().iter().count();
        assert_eq!(x, 2);
        for (id, key) in molecule.get_links().iter() {
            trace!(id,key.0, key.1);
        }
    }
  //  #[test]
//     fn test_parse_simple_kcf() {
//         let simple_kcf = r#"ATOM 4
// 1 1 1
// 1 1 2
// 1 1 2
// 1 1 2
// BONDS
// 1 2
// 2 3
// 3 4
// //"#;
//
//         let mut molecule = Molecule::new();
//         molecule.parse_simple_kcf(simple_kcf);
//
//         assert_eq!(molecule.atoms.len(), 4);
//         assert_eq!(molecule.links.len(), 3);
//
//         // Vérifier les atomes
//         assert_eq!(molecule.get_atom(1), Some((AtomType::C, AtomState::Normal)));
//         assert_eq!(molecule.get_atom(2), Some((AtomType::C, AtomState::Aromatic)));
//         assert_eq!(molecule.get_atom(3), Some((AtomType::C, AtomState::Aromatic)));
//         assert_eq!(molecule.get_atom(4), Some((AtomType::C, AtomState::Aromatic)));
//
//         // Vérifier les liens
//         assert_eq!(molecule.get_bond(1), Some((1, 2)));
//         assert_eq!(molecule.get_bond(2), Some((2, 3)));
//         assert_eq!(molecule.get_bond(3), Some((3, 4)));
//     }
//
//     #[test]
//     fn test_parse_simple_kcf_with_comments() {
//         let simple_kcf = r#"
// # Ceci est un commentaire
// ATOM 3
// 1 2 0  # Azote normal
// 2 1 1  # Carbone aromatique
// 3 3 0  # Oxygène normal
// BONDS
// 1 2
// 2 3
// # Fin des liens
// //
// "#;
//
//         let mut molecule = Molecule::new();
//         molecule.parse_simple_kcf(simple_kcf).unwrap();
//
//         assert_eq!(molecule.atoms.len(), 3);
//         assert_eq!(molecule.links.len(), 2);
//     }
//
//     #[test]
//     fn test_parse_simple_kcf_invalid() {
//         let invalid_kcf = r#"ATOM 2
// 1 1 1
// BONDS
// 1 3  # Atome 3 n'existe pas
// //"#;
//
//         let mut molecule = Molecule::new();
//         let result = molecule.parse_simple_kcf(invalid_kcf);
//         assert!(result.is_err());
//     }
//
//     #[test]
//     fn test_parse_kcf_auto_detection() {
//         // Format simple
//         let simple_content = r#"ATOM 2
// 1 1 1
// 2 1 1
// BONDS
// 1 2
// //"#;
//
//         let mut molecule = Molecule::new();
//         molecule.parse_kcf_auto(simple_content).unwrap();
//         assert_eq!(molecule.atoms.len(), 2);
//
//         // Format détaillé (votre format original)
//         let detailed_content = r#"ATOM 2
// 1 C8y C 1.0 2.0
// 2 O1a O 3.0 4.0
// BONDS
// 1 1 2 1
// //"#;
//
//         let mut molecule2 = Molecule::new();
//         molecule2.parse_kcf_auto(detailed_content).unwrap();
//         assert_eq!(molecule2.atoms.len(), 2);
//     }
}