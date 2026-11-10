#!/usr/bin/env python3
"""minimal_environment.py — build the smallest starting medium that lets
every reaction in a chemistry JSON fire at least once, ready to hand to
`simulate`.

This is the write side of the "scope" idea in reaction_network.py: rather
than asking what a *given* medium can produce, it asks for a minimal medium
that can produce (and thus fire) everything the reaction set is capable of.

Algorithm (greedy, not a global optimum — see reaction_network.py's docstring
for why an exact minimum generating set is a harder combinatorial problem):
  1. Repeatedly fire every reaction whose reactants are already available,
     closing under production.
  2. When nothing more can fire but some reaction is still short a reactant,
     seed exactly one missing reactant (from the reaction closest to firing)
     and resume. Every seed is therefore something no reaction in this file
     can produce on its own — it has to come from the starting medium.
  3. Repeat until every reaction's reactants are satisfiable.

A seed is either a free atom species ("a0") or a bonded pair ("a0.b1"). Free
atoms map to `-a KEY:COUNT`; bonded pairs can't be expressed with `-a` (there
is no CLI way to pre-link two atoms), so they are written out as a single
`.mol` file — one disjoint atom-pair fragment per bonded species, which
`simulate --mol` loads as one `Reactor::add_molecule` call.

Usage:
  python3 minimal_environment.py chem.json
  python3 minimal_environment.py chem.json --count 30 --bond-copies 10 --out-dir data/env
"""
import argparse
import os
import shlex

from ssa_chem import ENV, compute_scope, is_bond_species, load_reactions, symbol_to_key


def compute_minimal_seed(reactions):
    """Greedy minimal generating set: returns (seeds, reasons, scope).
    `seeds` is the ordered list of species that had to be supplied externally;
    `scope` is everything reachable (seeds + everything they cascade into)."""
    scope = {ENV}
    seeds = []
    reasons = {}

    def fire_pass():
        grew = False
        for r in reactions:
            if all(sp in scope for sp in r.reactants()):
                for sp in r.products():
                    if sp not in scope:
                        scope.add(sp)
                        grew = True
        return grew

    while True:
        while fire_pass():
            pass
        pending = [(r, {sp for sp in r.reactants() if sp not in scope}) for r in reactions]
        pending = [(r, missing) for r, missing in pending if missing]
        if not pending:
            break
        # Seed the reaction closest to firing first (fewest *distinct* missing
        # reactants — a self-reaction like a0+a0 must not count 'a0' twice).
        # On ties, prefer a reaction needing a free atom over one needing a
        # bond: atoms are the more generative unit (a single atom seed can
        # cascade via contact reactions into several different bonds), so
        # trying that path first keeps the final seed set smaller.
        r, missing = min(pending, key=lambda item: (len(item[1]),
                                                      all(is_bond_species(sp) for sp in item[1])))
        seed = next(iter(missing))
        scope.add(seed)
        seeds.append(seed)
        reasons[seed] = f"not producible by any reaction; needed by {r.to_string()} ({r.type_label()})"

    return seeds, reasons, scope


def build_mol_file(bond_seeds, copies):
    """One `.mol` block containing `copies` disjoint atom-pair fragments per
    bonded seed species — valid because Molecule/Reactor don't require the
    atoms in a .mol file to form a single connected component."""
    atom_lines, bond_lines = [], []
    next_id = 1
    for species in bond_seeds:
        a_sym, b_sym = is_bond_species(species)
        a_key, b_key = symbol_to_key(a_sym), symbol_to_key(b_sym)
        a_form, a_state = (a_key >> 8) & 0xFF, a_key & 0xFF
        b_form, b_state = (b_key >> 8) & 0xFF, b_key & 0xFF
        for _ in range(copies):
            aid, bid = next_id, next_id + 1
            next_id += 2
            atom_lines.append(f"{aid} {a_form} {a_state}")
            atom_lines.append(f"{bid} {b_form} {b_state}")
            bond_lines.append(f"{aid} {bid}")
    body = [f"ATOM {len(atom_lines)}"] + atom_lines + ["BONDS"] + bond_lines + ["//"]
    return "\n".join(body) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reactions_json", help="reactions JSON file")
    ap.add_argument("--count", type=int, default=20,
                     help="atom count per free-atom seed species (default 20)")
    ap.add_argument("--bond-copies", type=int, default=5,
                     help="pre-bonded pair copies per bonded seed species (default 5)")
    ap.add_argument("--out-dir", default="data/minimal_env",
                     help="where to write initial_bonds.mol and run_simulate.sh")
    ap.add_argument("--tmax", type=float, default=10000.0)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    reactions = load_reactions(args.reactions_json)
    seeds, reasons, scope = compute_minimal_seed(reactions)

    atom_seeds = [s for s in seeds if is_bond_species(s) is None]
    bond_seeds = [s for s in seeds if is_bond_species(s) is not None]

    print(f"{len(reactions)} reactions -> minimal seed set: {len(seeds)} species "
          f"({len(atom_seeds)} free atom(s), {len(bond_seeds)} pre-bonded pair(s))\n")
    for sp in seeds:
        kind = "bond" if is_bond_species(sp) else "atom"
        print(f"  [{kind}] {sp:<10} {reasons[sp]}")
    if not seeds:
        print("  (empty — every reaction is a production rule, or the file has no reactions)")

    # Sanity check: this seed set should make every reaction's reactants
    # satisfiable. compute_scope reruns the same closure independently.
    final_scope = compute_scope(reactions, seeds)
    stuck = [r for r in reactions if not all(sp in final_scope for sp in r.reactants())]
    if stuck:
        print(f"\nWARNING: {len(stuck)} reaction(s) still cannot fire from this seed set:")
        for r in stuck:
            print(f"  {r.to_string()}  (needs: {', '.join(r.reactants())})")
    else:
        print(f"\nAll {len(reactions)} reactions can fire from this seed set (verified).")

    os.makedirs(args.out_dir, exist_ok=True)

    atom_args = [f"0x{symbol_to_key(sp):04X}:{args.count}" for sp in atom_seeds]

    mol_path = None
    if bond_seeds:
        mol_path = os.path.join(args.out_dir, "initial_bonds.mol")
        with open(mol_path, "w") as f:
            f.write(build_mol_file(bond_seeds, args.bond_copies))
        print(f"\n{len(bond_seeds)} bonded species x {args.bond_copies} copies -> {mol_path}")

    cmd = ["cargo", "run", "--release", "--bin", "simulate", "--", "-r", args.reactions_json]
    for a in atom_args:
        cmd += ["-a", a]
    if mol_path:
        cmd += ["--mol", mol_path]
    cmd += ["--tmax", str(args.tmax), "--seed", str(args.seed),
            "-o", os.path.join(args.out_dir, "run")]

    cmd_str = " ".join(shlex.quote(c) for c in cmd)
    script_path = os.path.join(args.out_dir, "run_simulate.sh")
    with open(script_path, "w") as f:
        f.write(f"#!/bin/sh\n{cmd_str}\n")
    os.chmod(script_path, 0o755)

    print(f"\nReady-to-run command (also saved to {script_path}):")
    print(f"  {cmd_str}")


if __name__ == "__main__":
    main()
