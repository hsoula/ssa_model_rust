pub mod atoms;
pub mod reactions;
pub mod reactor;
pub mod links;
pub mod utils;
pub mod molecules;
mod event_log;
mod reactions_maps;
pub mod Gene8;
pub mod mol_genome;

pub use atoms::Atom;
pub use reactions::Reaction;
#[macro_export]
macro_rules! trace {
    ($($args: expr),*) => {
        print!("TRACE: file: {}, line: {}", file!(), line!());
        $(
            print!(", {}: {}", stringify!($args), $args);
        )*
        println!(); // Adds a final newline
    }
}
