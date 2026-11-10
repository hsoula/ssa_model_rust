#!/usr/bin/env python3
"""reaction_network.py — analyze a reactions JSON file produced by ssa_model.

Reads the same reactions-JSON format used by the Rust CLI tools
(gen_chemistry, evolve, Gene8::save_to_json, ...): a list of
{"s1","s2","p1","p2": "0xHHHH", "reaction_type": "0bBBBB", "reaction_speed": f64,
"tags": str?}, with the same wildcard convention as
`ssa_model::reactions::Reaction::news_from_json` (s1/s2 form byte == 0xFF means
"all forms").

Two things are computed, both purely from the rule set (no kinetics/SSA):

1. The reaction network graph: species nodes are either a free atom "a0"
   (form=a, state=0) or a bonded pair "a0.b1", plus one reaction node per
   rule (labeled like `Reaction::to_string()`, e.g. "a0+b1->a2.b0"), with
   edges reactant(s) -> reaction -> product(s).

2. The scope: given a starting medium (`--initial`), the set of species that
   can *possibly* be built by repeatedly firing rules whose reactants are
   already available, assuming an unbounded supply of anything already
   reached (the classic artificial-chemistry "scope"/closure — not a
   simulation, so reaction speed and stoichiometric counts are ignored).

Usage:
  python3 reaction_network.py chem.json
  python3 reaction_network.py chem.json --draw network.dot
  python3 reaction_network.py chem.json --initial 0x0100,0x0200 --draw-scope scope.png
"""
import argparse
import sys

import networkx as nx

from ssa_chem import ENV, compute_scope, load_reactions, parse_species_token

# ── Network graph ─────────────────────────────────────────────────────────

def build_graph(reactions):
    """Bipartite species/reaction graph: reactant(s) -> reaction -> product(s)."""
    graph = nx.MultiDiGraph()
    for idx, r in enumerate(reactions):
        rnode = f"R{idx}"
        graph.add_node(rnode, kind="reaction", label=r.to_string(),
                        rtype=r.type_label(), speed=r.speed)
        for sp in r.reactants():
            graph.add_node(sp, kind="species")
            graph.add_edge(sp, rnode)
        for sp in r.products():
            graph.add_node(sp, kind="species")
            graph.add_edge(rnode, sp)
    return graph


def species_nodes(graph):
    return sorted(n for n, d in graph.nodes(data=True) if d.get("kind") == "species")


def reaction_nodes(graph):
    return [n for n, d in graph.nodes(data=True) if d.get("kind") == "reaction"]


# ── Rendering ──────────────────────────────────────────────────────────────

def write_dot(graph, path, highlight=None):
    lines = ["digraph reaction_network {", "  rankdir=LR;", '  node [fontsize=10];']
    for n, data in graph.nodes(data=True):
        if data.get("kind") == "reaction":
            label = data.get("label", n).replace('"', r"\"")
            lines.append(f'  "{n}" [shape=box, style=filled, fillcolor=lightgoldenrod, label="{label}"];')
        else:
            if n == ENV:
                color = "lightgray"
            elif highlight is not None:
                color = "mediumseagreen" if n in highlight else "lightblue"
            else:
                color = "lightblue"
            lines.append(f'  "{n}" [shape=ellipse, style=filled, fillcolor={color}, label="{n}"];')
    for u, v in graph.edges():
        lines.append(f'  "{u}" -> "{v}";')
    lines.append("}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"DOT written to {path}")


def draw_matplotlib(graph, path, highlight=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pos = nx.spring_layout(graph, seed=42, k=0.9)
    species = species_nodes(graph)
    rxns = reaction_nodes(graph)

    fig, ax = plt.subplots(figsize=(max(8, len(species) * 0.6), max(6, len(species) * 0.4)))

    def species_color(n):
        if n == ENV:
            return "lightgray"
        if highlight is not None:
            return "mediumseagreen" if n in highlight else "lightblue"
        return "lightblue"

    nx.draw_networkx_nodes(graph, pos, nodelist=species, node_shape="o", node_size=600,
                            node_color=[species_color(n) for n in species], ax=ax)
    nx.draw_networkx_nodes(graph, pos, nodelist=rxns, node_shape="s", node_size=280,
                            node_color="lightgoldenrod", ax=ax)
    nx.draw_networkx_labels(graph, pos, labels={n: n for n in species}, font_size=8, ax=ax)
    nx.draw_networkx_edges(graph, pos, arrowsize=10, ax=ax,
                            connectionstyle="arc3,rad=0.05")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"Graph image written to {path}")


def write_graph(graph, path, highlight=None):
    if path.endswith(".dot"):
        write_dot(graph, path, highlight)
        print(f"  render with: dot -Tpng {path} -o {path[:-4]}.png")
        return
    try:
        draw_matplotlib(graph, path, highlight)
    except ImportError as e:
        fallback = path.rsplit(".", 1)[0] + ".dot"
        print(f"matplotlib unavailable ({e}); writing {fallback} instead", file=sys.stderr)
        write_dot(graph, fallback, highlight)
        print(f"  render with: dot -Tpng {fallback} -o {path}")


# ── CLI ──────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reactions_json", help="reactions JSON file")
    ap.add_argument("--initial",
                     help="comma-separated starting medium: hex atom keys (0x0100) "
                          "and/or symbolic species (a0, a0.b1); enables scope computation")
    ap.add_argument("--draw", metavar="FILE",
                     help="write the full reaction network (.dot always works; "
                          ".png/.pdf/.svg via matplotlib if installed)")
    ap.add_argument("--draw-scope", metavar="FILE",
                     help="like --draw, restricted to the species/reactions reachable "
                          "in the computed scope (requires --initial)")
    args = ap.parse_args()

    reactions = load_reactions(args.reactions_json)
    graph = build_graph(reactions)
    species = species_nodes(graph)

    print(f"Loaded {len(reactions)} reactions -> {len(species)} species "
          f"({len(reaction_nodes(graph))} reaction nodes)\n")

    print("Reactions:")
    for i, r in enumerate(reactions):
        print(f"  [{i:3}] {r.to_string():<16} {r.type_label():<12} speed={r.speed:.3f}")

    print(f"\nSpecies ({len(species)}): {', '.join(species)}")

    scope = None
    if args.initial:
        initial = [parse_species_token(t) for t in args.initial.split(",") if t.strip()]
        scope = compute_scope(reactions, initial)
        unreachable = [s for s in species if s not in scope]
        print(f"\n=== Scope from initial medium {{{', '.join(initial)}}} ===")
        print(f"Reachable  ({len(scope)}): {', '.join(sorted(scope))}")
        print(f"Unreachable ({len(unreachable)}): {', '.join(unreachable) if unreachable else '(none)'}")

    if args.draw:
        write_graph(graph, args.draw)

    if args.draw_scope:
        if scope is None:
            sys.exit("--draw-scope requires --initial")
        active_reactions = {
            n for n in reaction_nodes(graph)
            if all(pred in scope for pred in graph.predecessors(n))
        }
        keep = set(scope) | active_reactions | {ENV}
        sub = graph.subgraph(keep)
        write_graph(sub, args.draw_scope, highlight=scope)


if __name__ == "__main__":
    main()
