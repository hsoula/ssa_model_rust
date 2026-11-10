#!/usr/bin/env python3
"""compound_reaction_graph.py — visualize KEGG compound graphs and reactions.

Two local datasets, independent of the Rust ssa_model chemistry:

  molecules/json/<CID>.json   — one KCF-style atom/bond graph per KEGG compound
                                 (entry_id, formula, atoms:[{index,type,symbol,x,y}],
                                 bonds:[{atom1,atom2,bond_type}])
  reactions/all_reactions.json — KEGG reaction records (entry, name, equation,
                                 substrates:[CID...], products:[CID...])

This script builds an actual atom-bond graph (networkx) per compound — using
each atom's own (x,y) from the KCF record, so the layout looks like a real
skeletal structure rather than a generic force-directed blob — and, for a
reaction, tiles the substrate graphs on the left and product graphs on the
right of an arrow. It also reports basic graph-theoretic properties per
compound (atom/bond/ring counts, element composition) and a best-effort
element-balance check across a reaction, computed from the graphs themselves
rather than the textual formula — a first step towards inferring chemistry
properties from graph structure rather than hand-curated chemical rules.

Usage:
  python3 compound_reaction_graph.py --reaction R00001
  python3 compound_reaction_graph.py --reaction 1 --draw R00001.dot
  python3 compound_reaction_graph.py --compound C00022 --draw C00022.png
  python3 compound_reaction_graph.py --all
  python3 compound_reaction_graph.py --all --export-dir out/reactions --limit 50
"""
import argparse
import json
import os
import re
import sys
from collections import Counter

import networkx as nx

# ── Appearance ───────────────────────────────────────────────────────────────

ELEMENT_COLORS = {
    "C": "#4d4d4d", "O": "#e0301e", "N": "#3050f8", "P": "#ff8000",
    "S": "#ffd123", "H": "#dddddd", "Cl": "#1ff01f", "F": "#90e050",
    "Br": "#a62929", "I": "#940094", "Na": "#ab5cf2", "Mg": "#8aff00",
    "K": "#8f40d4", "Fe": "#e06633", "Zn": "#7d80b0",
}
DEFAULT_COLOR = "#c39bd3"  # generic/unknown KCF placeholders (R, *, R#, ...)
BOND_WIDTH = {"1": 1.4, "2": 2.6, "3": 3.8}

# KCF coordinates use ~1.4 units per bond length. Graphviz's pinned "x,y!"
# positions are interpreted in points (72/inch), so raw KCF units would place
# atoms a fraction of a point apart; scale them up to a sane on-page spacing.
LAYOUT_SCALE = 36.0


# ── Compound graphs ──────────────────────────────────────────────────────────

MAX_COMPOUND_FILE_BYTES = 2_000_000  # typical compound JSON is a few KB; guards against corrupt files


def load_compound(compound_id, mol_dir):
    path = os.path.join(mol_dir, f"{compound_id}.json")
    if not os.path.exists(path):
        return None
    try:
        if os.path.getsize(path) > MAX_COMPOUND_FILE_BYTES:
            print(f"  warning: skipping {path} (unexpectedly large, likely corrupt)", file=sys.stderr)
            return None
        with open(path, encoding="utf-8", errors="replace") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        print(f"  warning: could not read {path}: {e}", file=sys.stderr)
        return None


def compound_graph(mol_json):
    """KCF-style {atoms, bonds} JSON -> nx.Graph with element/symbol/x/y and
    bond_type node/edge attributes."""
    g = nx.Graph(name=mol_json.get("entry_id", "?"))
    for a in mol_json.get("atoms", []):
        g.add_node(a["index"], element=a.get("type") or a.get("symbol", "?"),
                   symbol=a.get("symbol", ""),
                   x=a.get("x", 0.0) * LAYOUT_SCALE, y=a.get("y", 0.0) * LAYOUT_SCALE)
    for b in mol_json.get("bonds", []):
        g.add_edge(b["atom1"], b["atom2"], bond_type=b.get("bond_type", "1"))
    return g


def graph_stats(g):
    n, m = g.number_of_nodes(), g.number_of_edges()
    comps = nx.number_connected_components(g) if n else 0
    rings = max(m - n + comps, 0) if n else 0  # cyclomatic number (first Betti number)
    degrees = [d for _, d in g.degree()]
    return {
        "n_atoms": n,
        "n_bonds": m,
        "n_rings": rings,
        "n_components": comps,
        "avg_degree": (sum(degrees) / n) if n else 0.0,
        "elements": Counter(nx.get_node_attributes(g, "element").values()),
    }


# ── Reactions ────────────────────────────────────────────────────────────────

def load_reactions(path):
    with open(path) as f:
        raw = json.load(f)
    return {r["entry"]: r for r in raw}


def normalize_reaction_id(query, reactions):
    q = query.strip().upper()
    if q in reactions:
        return q
    if q.isdigit() and f"R{int(q):05d}" in reactions:
        return f"R{int(q):05d}"
    raise KeyError(f"reaction '{query}' not found")


def parse_side_coeffs(side_text):
    """'16 C00002 + 16 C00001' -> {'C00002': 16, 'C00001': 16}.
    KEGG separates top-level terms with ' + ' (spaces), so this does not
    confuse the terms with the algebraic '+' inside e.g. '(n+1)'. Symbolic
    coefficients ('n', '(n+1)', '(n-1)', ...) map to None (variable)."""
    out = {}
    for term in side_text.split(" + "):
        m = re.match(r"^(.*?)\s*(C\d{5})$", term.strip())
        if not m:
            continue
        coeff_str, cid = m.group(1).strip(), m.group(2)
        if coeff_str == "":
            coeff = 1
        elif coeff_str.isdigit():
            coeff = int(coeff_str)
        else:
            coeff = None
        out[cid] = coeff
    return out


def parse_equation(equation):
    sep = "<=>" if "<=>" in equation else "=>"
    left, _, right = equation.partition(sep)
    return parse_side_coeffs(left), parse_side_coeffs(right)


def element_balance(compound_coeffs, mol_dir):
    """Sum(coeff * element_count) across a reaction side, from the graphs
    themselves. Returns (totals, symbolic, complete) — `complete` is False if
    any compound's molecule JSON is missing locally."""
    total, symbolic = Counter(), False
    for cid, coeff in compound_coeffs.items():
        mol = load_compound(cid, mol_dir)
        if mol is None:
            return total, symbolic, False
        if coeff is None:
            symbolic, coeff = True, 1
        for el, cnt in Counter(nx.get_node_attributes(compound_graph(mol), "element").values()).items():
            total[el] += cnt * coeff
    return total, symbolic, True


# ── Layout: tile each compound's own (x,y) atom coordinates side by side ────

def layout_side_by_side(substrate_graphs, product_graphs, gap=3.0 * LAYOUT_SCALE):
    """Returns (combined_graph, pos, centers, reaction_x). Node ids in the
    combined graph are (compound_id, atom_index) to keep compounds disjoint."""
    combined = nx.Graph()
    pos, centers = {}, {}
    cursor = 0.0

    def place(graphs):
        nonlocal cursor
        for cid, g in graphs:
            if g.number_of_nodes() == 0:
                continue
            xs = [d["x"] for _, d in g.nodes(data=True)]
            ys = [d["y"] for _, d in g.nodes(data=True)]
            x0, y0 = min(xs), sum(ys) / len(ys)
            width = max(xs) - x0
            for n, d in g.nodes(data=True):
                key = (cid, n)
                combined.add_node(key, **d)
                pos[key] = (cursor + (d["x"] - x0), d["y"] - y0)
            for u, v, d in g.edges(data=True):
                combined.add_edge((cid, u), (cid, v), **d)
            centers[cid] = cursor + width / 2
            cursor += width + gap

    place(substrate_graphs)
    reaction_x = cursor + gap / 2 if product_graphs else None
    if product_graphs:
        cursor += gap * 2
    place(product_graphs)
    return combined, pos, centers, reaction_x


def _load_side(compound_ids, mol_dir, missing):
    graphs = []
    for cid in compound_ids:
        mol = load_compound(cid, mol_dir)
        if mol is None:
            missing.append(cid)
            continue
        graphs.append((cid, compound_graph(mol)))
    return graphs


# ── Rendering ────────────────────────────────────────────────────────────────

def _draw_matplotlib(combined, pos, centers, reaction_x, title, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    width = max(8.0, 2.5 * max(len(centers), 1))
    fig, ax = plt.subplots(figsize=(width, 5))

    widths = [BOND_WIDTH.get(d.get("bond_type", "1"), 1.4) for _, _, d in combined.edges(data=True)]
    nx.draw_networkx_edges(combined, pos, ax=ax, width=widths or 1.0)

    colors = [ELEMENT_COLORS.get(d["element"], DEFAULT_COLOR) for _, d in combined.nodes(data=True)]
    nx.draw_networkx_nodes(combined, pos, node_color=colors, node_size=140,
                            edgecolors="black", linewidths=0.4, ax=ax)

    labels = {n: d["element"] for n, d in combined.nodes(data=True) if d["element"] != "C"}
    nx.draw_networkx_labels(combined, pos, labels=labels, font_size=7, ax=ax)

    for cid, cx in centers.items():
        ys = [pos[n][1] for n in combined.nodes if n[0] == cid]
        ax.text(cx, max(ys, default=0) + LAYOUT_SCALE, cid, ha="center", fontsize=9, fontweight="bold")

    if reaction_x is not None:
        ax.annotate("", xy=(reaction_x + LAYOUT_SCALE / 3, 0), xytext=(reaction_x - LAYOUT_SCALE / 3, 0),
                    arrowprops=dict(arrowstyle="-|>", lw=2))

    ax.set_title(title, fontsize=10)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def _write_dot(combined, pos, centers, reaction_x, title, path):
    lines = ["graph reaction {", "  layout=neato;",
             '  node [shape=circle, fontsize=8, style=filled, fixedsize=true, width=0.22];']
    for n, d in combined.nodes(data=True):
        x, y = pos[n]
        nid = f'"{n[0]}::{n[1]}"'
        color = ELEMENT_COLORS.get(d["element"], DEFAULT_COLOR)
        lines.append(f'  {nid} [pos="{x:.2f},{y:.2f}!", label="{d["element"]}", fillcolor="{color}"];')
    for u, v, d in combined.edges(data=True):
        w = {"1": 1, "2": 2, "3": 3}.get(d.get("bond_type", "1"), 1)
        lines.append(f'  "{u[0]}::{u[1]}" -- "{v[0]}::{v[1]}" [penwidth={w}];')

    # Compound-id label above each molecule cluster
    for cid, cx in centers.items():
        ys = [pos[n][1] for n in combined.nodes if n[0] == cid]
        top = max(ys, default=0) + LAYOUT_SCALE
        lines.append(f'  "label::{cid}" [pos="{cx:.2f},{top:.2f}!", shape=plaintext, '
                     f'fixedsize=false, label="{cid}", fontsize=12];')

    # Reactants -> products arrow
    if reaction_x is not None:
        lines.append(f'  "arrow" [pos="{reaction_x:.2f},0!", shape=plaintext, '
                     f'fixedsize=false, label="→", fontsize=28];')

    safe_title = title.replace('"', "'")
    lines.append(f'  label="{safe_title}"; labelloc=t; fontsize=14;')
    lines.append("}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_diagram(combined, pos, centers, reaction_x, title, path, warn=print):
    if path.endswith(".dot"):
        _write_dot(combined, pos, centers, reaction_x, title, path)
        warn(f"DOT written to {path}  (render with: neato -n -Tpng {path} -o {path[:-4]}.png)")
        return
    try:
        _draw_matplotlib(combined, pos, centers, reaction_x, title, path)
        warn(f"Image written to {path}")
    except ImportError as e:
        fallback = path.rsplit(".", 1)[0] + ".dot"
        warn(f"matplotlib unavailable ({e}); writing {fallback} instead")
        _write_dot(combined, pos, centers, reaction_x, title, fallback)
        warn(f"  render with: neato -n -Tpng {fallback} -o {path}")


# ── High-level actions ───────────────────────────────────────────────────────

def analyze_compound(cid, mol_dir, draw_path):
    mol = load_compound(cid, mol_dir)
    if mol is None:
        sys.exit(f"no local molecule JSON for {cid} in {mol_dir}")
    g = compound_graph(mol)
    st = graph_stats(g)
    name = (mol.get("names") or [""])[0]
    print(f"{cid}  {mol.get('formula', '')}  {name}")
    print(f"  atoms={st['n_atoms']}  bonds={st['n_bonds']}  rings={st['n_rings']}  "
          f"components={st['n_components']}  avg_degree={st['avg_degree']:.2f}")
    print(f"  elements: {dict(st['elements'])}")
    if draw_path:
        combined, pos, centers, reaction_x = layout_side_by_side([(cid, g)], [])
        write_diagram(combined, pos, centers, reaction_x, f"{cid}  {name}", draw_path)


def analyze_reaction(reaction, mol_dir, draw_path):
    subs, prods = reaction.get("substrates", []), reaction.get("products", [])
    print(f"{reaction['entry']}  {reaction.get('name', '')}")
    print(f"  equation: {reaction.get('equation', '')}")
    print(f"  {' + '.join(subs)}  ->  {' + '.join(prods)}")

    print("\n  Per-compound graph stats:")
    for cid in dict.fromkeys(subs + prods):
        mol = load_compound(cid, mol_dir)
        if mol is None:
            print(f"    {cid:<10} (no local molecule JSON)")
            continue
        st = graph_stats(compound_graph(mol))
        print(f"    {cid:<10} {mol.get('formula', ''):<14} atoms={st['n_atoms']:<3} "
              f"bonds={st['n_bonds']:<3} rings={st['n_rings']}")

    left_c, right_c = parse_equation(reaction.get("equation", ""))
    left_sum, left_sym, left_ok = element_balance(left_c or {cid: 1 for cid in subs}, mol_dir)
    right_sum, right_sym, right_ok = element_balance(right_c or {cid: 1 for cid in prods}, mol_dir)
    if left_ok and right_ok:
        status = "balanced" if left_sum == right_sum else "UNBALANCED"
        note = " (symbolic stoichiometry present, coeff=1 assumed)" if (left_sym or right_sym) else ""
        print(f"\n  Element balance from graphs (heavy atoms only — KCF omits implicit H): {status}{note}")
        print(f"    substrates: {dict(left_sum)}")
        print(f"    products:   {dict(right_sum)}")
    else:
        print("\n  Element balance: skipped (missing local molecule JSON for one or more compounds)")

    if draw_path:
        missing = []
        sub_graphs = _load_side(subs, mol_dir, missing)
        prod_graphs = _load_side(prods, mol_dir, missing)
        if missing:
            print(f"  (no local molecule JSON for: {', '.join(missing)}; drawing the rest)")
        if not sub_graphs and not prod_graphs:
            print("  nothing to draw"); return
        combined, pos, centers, reaction_x = layout_side_by_side(sub_graphs, prod_graphs)
        title = f"{reaction['entry']} — {reaction.get('name', '')}"
        write_diagram(combined, pos, centers, reaction_x, title, draw_path)


def summarize_all(reactions, mol_dir, export_dir, limit, fmt):
    covered, shapes = [], Counter()
    for r in reactions.values():
        subs, prods = r.get("substrates", []), r.get("products", [])
        ok = all(os.path.exists(os.path.join(mol_dir, f"{c}.json")) for c in subs + prods)
        if ok:
            covered.append(r)
        shapes[(len(subs), len(prods))] += 1

    n = len(reactions)
    print(f"{n} reactions total; {len(covered)} fully covered by local molecule JSONs "
          f"({100 * len(covered) / n:.1f}%)")
    print("Most common (#substrates, #products) shapes:")
    for shape, cnt in shapes.most_common(10):
        print(f"  {shape}: {cnt}")

    if export_dir:
        os.makedirs(export_dir, exist_ok=True)
        todo = covered if limit == 0 else covered[:limit]
        print(f"\nExporting {len(todo)} reaction diagram(s) to {export_dir}/ (.{fmt}) ...")
        for r in todo:
            out = os.path.join(export_dir, f"{r['entry']}.{fmt}")
            try:
                missing = []
                sub_graphs = _load_side(r.get("substrates", []), mol_dir, missing)
                prod_graphs = _load_side(r.get("products", []), mol_dir, missing)
                combined, pos, centers, reaction_x = layout_side_by_side(sub_graphs, prod_graphs)
                title = f"{r['entry']} — {r.get('name', '')}"
                write_diagram(combined, pos, centers, reaction_x, title, out, warn=lambda *a: None)
            except Exception as e:
                print(f"  {r['entry']}: failed ({e})")
        print("Done.")


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reactions-json", default="reactions/all_reactions.json")
    ap.add_argument("--molecules-dir", default="molecules/json")
    ap.add_argument("--reaction", help="reaction id (R00001) or bare number (1)")
    ap.add_argument("--compound", help="a single compound id (e.g. C00022), no reaction context")
    ap.add_argument("--all", action="store_true", help="scan every reaction: coverage + shape stats")
    ap.add_argument("--draw", metavar="FILE",
                     help="diagram output for --reaction/--compound "
                          "(.dot always works; .png/.svg/... via matplotlib if installed)")
    ap.add_argument("--export-dir", metavar="DIR",
                     help="with --all: batch-render every fully-covered reaction into DIR")
    ap.add_argument("--limit", type=int, default=20,
                     help="cap on reactions rendered by --export-dir (0 = no limit, default 20)")
    ap.add_argument("--format", default="dot", choices=["dot", "png", "svg", "pdf"],
                     help="file format used by --export-dir (default: dot)")
    args = ap.parse_args()

    if not (args.reaction or args.compound or args.all):
        ap.error("specify one of --reaction, --compound, --all")

    if args.compound:
        analyze_compound(args.compound, args.molecules_dir, args.draw)
        return

    reactions = load_reactions(args.reactions_json)

    if args.reaction:
        try:
            rid = normalize_reaction_id(args.reaction, reactions)
        except KeyError as e:
            sys.exit(str(e))
        analyze_reaction(reactions[rid], args.molecules_dir, args.draw)
        return

    summarize_all(reactions, args.molecules_dir, args.export_dir, args.limit, args.format)


if __name__ == "__main__":
    main()
