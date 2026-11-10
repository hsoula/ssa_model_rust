# ssa_model

A stochastic simulation algorithm (Gillespie SSA) for an artificial chemistry in
the style of Tim Hutton's molecule/atom experiments: a small alphabet of typed
"atoms" that link into molecule graphs under a set of probabilistic reaction
rules. On top of the simulator sit two genetic-algorithm encodings (`Gene8` for
reaction chemistries, `MolGenome` for molecules) used to *evolve* chemistries
rather than hand-write them.

## Core model

Everything is built from five pieces in `src/`:

| Module | Role |
|---|---|
| [`atoms.rs`](src/atoms.rs) | `Atom`: an id + `form` (type, 0–255) + `state` (0–255) + a set of link ids to other atoms. |
| [`links.rs`](src/links.rs) | `Link`: an undirected bond between two atom ids, keyed by the pair's `(form,state)` values for fast lookup. |
| [`molecules.rs`](src/molecules.rs) | `Molecule`: a standalone atoms+bonds graph (petgraph-backed). Load/save `.mol` text and KEGG-style JSON, generate random molecules, and compare molecules (isomorphism, composition/degree/Weisfeiler-Lehman distances). |
| [`reactions.rs`](src/reactions.rs) | `Reaction`: a rule `(substrate_a, substrate_b) → (product_a, product_b)` plus a `reaction_type` bitflag and a speed constant. Loads from JSON, including "wildcard" reactions that expand into one rule per form/state. |
| [`reactor.rs`](src/reactor.rs) | `Reactor`: the SSA engine itself — holds the live atom/link population and reaction list, computes propensities, and runs the Gillespie main loop. |

Two more modules add genome encodings on top of this model:

| Module | Role |
|---|---|
| [`Gene8.rs`](src/Gene8.rs) | Linear hex-nibble genome that decodes to a list of `Reaction`s. Used by the GA binaries to evolve a *chemistry*. Full token table in [`docs/gene8_encoding.md`](docs/gene8_encoding.md) (regenerate with `gene_doc`). |
| [`mol_genome.rs`](src/mol_genome.rs) | Linear hex-nibble genome that decodes to a list of `Molecule`s (linear atom chains). Used to evolve/seed *molecules*. Encoding rules are documented in the module's doc comment. |

Supporting modules: [`utils.rs`](src/utils.rs) (`AtomType`/`AtomState` enums, KCF-symbol parsing, the global `ObjectId` generator, `JsonReaction` (de)serialization), [`event_log.rs`](src/event_log.rs) (counts tagged reaction firings by name) and [`reactions_maps.rs`](src/reactions_maps.rs) (name ⇄ id registry for tags).

### Key encoding: the atom "key"

An atom's `(form, state)` pair is packed into a `u16` key = `form << 8 | state`,
and printed as a letter+digit symbol (`form + 'a'-1`, e.g. form 1/state 2 → `a2`).
A `Reaction`'s substrate/product pair is packed the same way into a `u32`
(`key_a << 16 | key_b`). This is the encoding used throughout the reactor,
`Gene8`, and every CLI tool's console output.

### Atom specs (`-a KEY:COUNT`) and wildcards

`simulate` and `evolve` take initial atom populations as `-a KEY:COUNT`, repeatable.
`KEY` is the 4-hex-digit `0xFFSS` byte pair above. Either byte's 2-digit group
may contain a `*` to wildcard it (parsed by `ssa_model::utils::parse_atom_spec`):

- `0x0*0*:50` — wildcard both bytes → every form (`1..=max_form`) × every state (`0..=max_state`), 50 atoms each.
- `0x0*03:50` — wildcard the form byte only → every form, state fixed at `0x03`.
- `0x01**:50` — wildcard the state byte only → form fixed at `0x01`, every state.

`max_form`/`max_state` default to 10 and are set with `--max-form`/`--max-state`.
Plain concrete keys (`0x0100:300`, decimal `256:300`) are unaffected.

### Reaction types

`reaction_type` is a 4-bit flag (`ReactionFlag` in `reactions.rs`):

- **Production** (`0b0100`): a bare atom is spawned (∅ → atom), propensity = speed.
- **Degradation** (`0b1000`): a matching atom is removed (atom → ∅), propensity = speed × count.
- **Substrate-contact** (bit `0b10` set): two *free* (unlinked) atoms collide; the reaction either bonds them (`product_contact = 0`) or just changes their states (`product_contact = 1`).
- **Link reaction** (substrate-contact bit clear): acts on an *existing* bond; either breaks it (`product_contact = 1`) or just changes both atoms' states.

### The Gillespie loop (`Reactor::main_loop`)

Standard direct-method SSA: compute a propensity per reaction
(`get_propensity_vector`), draw the time to the next event from the total
propensity `a0`, pick which reaction fires proportionally to its share of `a0`,
apply it (`apply_reaction` — spawn/remove an atom, bond/unbond a pair, or
recolor states), then repeat until `tmax` or `a0 ≈ 0`. The reactor also keeps
a fire-count histogram, a full `(time, reaction_idx)` sequence, and a
snapshot of atom/link counts after every step — enough to reconstruct a run
for a report (`to_markdown` / `save_markdown`) without re-simulating.

## Binaries

All binaries are `clap`-based CLIs (run with `cargo run --release --bin <name> -- --help`
for the full flag list) except `simple_ga` and `gene_doc`, which are fixed
single-purpose scripts. They fall into four groups:

**Build a chemistry / molecule set**
- `gen_chemistry` — generate a random reaction set as JSON (+ a readable `.md` report alongside it).
- `gen_molecules` — generate random connected molecules (Prüfer-tree + optional extra edges) as JSON + `.mol`.

**Run a simulation**
- `simulate` — load reaction JSON file(s) and an initial population (`-a KEY:COUNT` atoms and/or a `--mol` file of pre-bonded pairs — at least one of the two is required), run the SSA to `tmax`, dump surviving molecules and a `summary.json`.

**Evolve a chemistry with a GA** (population of `Gene8` genomes, tournament selection, crossover + point mutation)
- `evolve` — generic fitness (rewards linked-atom count, size diversity, max molecule size).
- `evolve_target` — fitness scored against one or more target molecules (composition similarity, size proximity, exact isomorphism), remapping `Gene8` forms onto the targets' actual atom forms.
- `simple_ga` — a fixed single-form/three-state micro-experiment (target = 5-atom chain) that additionally logs full phylogeny (parent ids, point-mutation positions) and per-reaction state snapshots, meant to be paired with `report_run`.

**Inspect / convert / report**
- `genome_test` — decode a random or fixed-length `Gene8` genome and print/save its valid reactions.
- `gene_doc` — print the `Gene8` token table + worked examples (regenerates `docs/gene8_encoding.md`).
- `mol_genome_test` — decode a `MolGenome` (hex string or random) and print the molecules it encodes.
- `mol_genome_reactor` — decode a `MolGenome`, spawn N copies of each molecule plus free atoms into a `Reactor`, load an optional chemistry JSON, and run the simulation.
- `mol_distances` — load all molecule JSONs from a directory and write pairwise composition / degree-sequence / Weisfeiler-Lehman distance matrices as CSV.
- `reactions_to_md` — render a reactions JSON file as a markdown table.
- `report_run` — re-simulate one individual from a `simple_ga` output directory (by id) and write a self-contained markdown report.

### Typical pipeline

```bash
# 1. Generate a random chemistry
cargo run --release --bin gen_chemistry -- -n 20 -f 4 -s 4 --seed 42 -o data/my_chem.json

# 2. Run it
cargo run --release --bin simulate -- \
  -r data/my_chem.json \
  -a 0x0100:300 -a 0x0200:300 -a 0x0300:300 \
  --tmax 10000 --seed 42 -o data/output/run1

# 3. Or evolve a chemistry instead of hand-generating one
cargo run --release --bin evolve -- \
  -a 0x0100:200 -a 0x0200:200 \
  --population 30 --generations 100 --tmax 500 --seed 42 -o data/evolve/run1

# 4. Evolve toward a specific target molecule
cargo run --release --bin evolve_target -- \
  -m molecules/json/C01101.json \
  --n-atoms 100 --genome-size 300 --population 30 --generations 200 --tmax 2000 \
  --seed 42 -o data/evolve/L-Ribulose
```

More end-to-end examples are in [`scripts.sh`](scripts.sh) and
[`generation_simulation.sh`](generation_simulation.sh).

## Data formats

- **`.mol`** — plain-text molecule format: `ATOM <n>` then `<id> <form> <state>` per line, `BONDS` then `<atom_a> <atom_b>` per line, terminated by `//`. One or more molecules can be concatenated in a single file.
- **Molecule JSON** — `{"name", "n_atoms", "n_bonds", "atoms": [{"id","form","state"}...], "bonds": [{"id","a","b"}...]}`, or the KEGG-style `{"entry_id","formula","atoms":[{"index","type","symbol"}...],"bonds":[{"index","atom1","atom2"}...]}` accepted by `Molecule::from_json_value`.
- **Reactions JSON** — array of `{"s1","s2","p1","p2": "0xHHHH", "reaction_type": "0bBBBB", "reaction_speed": f64, "tags": string?}`. A substrate field of `0xFF00` (or `0xFE00` on `s2` only) marks a *wildcard* reaction, expanded by `Reaction::news_from_json` into one concrete rule per form (and, for `0xFF00`/`0xFF00`, per form × state).

## Reaction network & scope (Python)

[`tools/reaction_network.py`](tools/reaction_network.py) analyzes a reactions
JSON file directly, without running any simulation:

- **Network graph** — a bipartite graph of species (`a0` free atoms, `a0.b1`
  bonds) and reactions (`a0+b1->a2.b0`-style nodes), written as Graphviz
  `.dot` (always works) or `.png`/`.svg`/... via `matplotlib` if installed.
- **Scope** — given a starting medium (`--initial`), the fixed-point closure
  of species reachable by repeatedly firing any reaction whose reactants are
  already present, assuming an unbounded supply of anything already reached.
  This is a pure reachability computation (no kinetics, no SSA, no counts) —
  it answers "what could possibly appear", not "what will appear".

```bash
# requires: pip install networkx  (matplotlib optional — falls back to .dot)
python3 tools/reaction_network.py data/my_chem.json
python3 tools/reaction_network.py data/my_chem.json --draw network.dot
python3 tools/reaction_network.py data/my_chem.json \
  --initial 0x0100,0x0200 --draw-scope scope.dot
dot -Tpng network.dot -o network.png   # render a .dot with Graphviz
```

[`tools/minimal_environment.py`](tools/minimal_environment.py) is the write
side of the same closure computation: instead of checking what a *given*
medium can reach, it derives the smallest starting medium that lets *every*
reaction in the file fire at least once. It greedily seeds only species no
reaction can produce (closing under production after each seed), so anything
reachable by cascading through the given reactions is never seeded directly.
Seeds are either free atoms (`-a KEY:COUNT`) or bonded pairs — since `-a` has
no way to pre-link two atoms, bonded seeds are written to a single `.mol`
file (multiple disjoint atom-pair fragments in one file, which `simulate
--mol` loads in one shot). It emits a ready-to-run `simulate` command:

```bash
python3 tools/minimal_environment.py data/my_chem.json --out-dir data/env
sh data/env/run_simulate.sh
```

(`tools/ssa_chem.py` holds the reaction-JSON parsing shared by both scripts.)

## Real KEGG compound/reaction graphs (Python)

Independent of the abstract atom/state chemistry above, [`molecules/json/`](molecules/json)
holds one real KEGG compound record per file (`entry_id`, `formula`,
`atoms:[{index,type,symbol,x,y}]`, `bonds:[{atom1,atom2,bond_type}]` — the
same KCF-style shape `Molecule::from_json_value` parses in the Rust code),
and [`reactions/all_reactions.json`](reactions/all_reactions.json) holds real
KEGG reaction records (`entry`, `name`, `equation`, `substrates`, `products`
as compound ids). [`tools/compound_reaction_graph.py`](tools/compound_reaction_graph.py)
turns these into actual molecule/reaction diagrams and basic graph-theoretic
properties — a first step towards inferring chemistry from graph structure
rather than curated rules:

```bash
# requires: pip install networkx  (matplotlib optional — falls back to Graphviz .dot)
python3 tools/compound_reaction_graph.py --reaction R00006 --draw r00006.dot
neato -n -Tpng r00006.dot -o r00006.png     # render a .dot with Graphviz

python3 tools/compound_reaction_graph.py --compound C00022     # one molecule's graph + stats
python3 tools/compound_reaction_graph.py --all                 # coverage + shape stats over all reactions
python3 tools/compound_reaction_graph.py --all --export-dir out/ --limit 50
```

For a reaction it tiles each substrate/product's own atom graph (using the
compound's real 2D atom coordinates, so it reads like an actual skeletal
structure) left/right of an arrow, and reports per-compound atom/bond/ring
counts plus a best-effort *heavy-atom* element balance computed straight from
the graphs (KCF omits implicit hydrogens, so this is not a full mass balance;
symbolic stoichiometry like KEGG's `n`/`(n+1)` polymer notation is flagged
and treated as coefficient 1). Note: `molecules/json/C00001.json` (water) is
corrupt in this checkout (171MB of mostly null bytes instead of a small JSON
record) — the script detects and skips oversized/unparseable compound files
rather than failing, but any reaction balance/diagram involving C00001 will
be missing that side.

## Tests

`tests/*.rs` cover `Atom`, `Molecule`, `Reaction`/`Reactor` behaviour and the
`MolGenome` decoder has its own unit tests inline in `src/mol_genome.rs`. Run
with `cargo test`.

## Repository layout

```
src/            core library (lib.rs re-exports the modules above)
src/bin/        the 13 CLI binaries described above
tools/          Python analysis scripts (ssa_chem.py, reaction_network.py,
                minimal_environment.py, compound_reaction_graph.py)
tests/          integration tests for the core library
chemistries/    example generated chemistry JSON files
molecules/json/ real KEGG compound graphs (also used as evolve_target targets)
reactions/      real KEGG reaction records (all_reactions.json)
data/, results/ simulation outputs (generated, gitignored)
docs/           generated reference docs (gene8_encoding.md)
```
