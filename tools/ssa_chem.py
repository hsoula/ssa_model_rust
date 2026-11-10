"""ssa_chem.py — shared reaction-JSON parsing for the ssa_model Python tools.

Species/reaction semantics mirror src/reactions.rs and src/reactor.rs:
a species is either a free atom "a0" (form=a, state=0) or a bonded pair
"a0.b1"; a reaction's reactants/products are derived from its reaction_type
flag exactly as Reactor::apply_reaction interprets it.
"""
import json

ENV = "∅"  # environment: source of production, sink of degradation

# ── Reaction type flags (mirrors ssa_model::reactions::ReactionFlag) ────────
SUBSTRATE_CONTACT = 0b0010
PRODUCT_CONTACT   = 0b0001
PRODUCTION        = 0b0100
DEGRADATION       = 0b1000


def form_symbol(form: int) -> str:
    return chr((form & 0xFF) + 96)


def atom_symbol(key: int) -> str:
    """u16 atom key (form<<8 | state) -> display symbol, e.g. 'a2'."""
    form, state = (key >> 8) & 0xFF, key & 0xFF
    return f"{form_symbol(form)}{state}"


def bond_symbol(key_a: int, key_b: int) -> str:
    """Undirected bond species, canonicalised so a0.b1 == b1.a0."""
    return ".".join(sorted((atom_symbol(key_a), atom_symbol(key_b))))


def symbol_to_key(symbol: str) -> int:
    """Inverse of atom_symbol: 'a2' -> 0x0102. Raises ValueError if malformed."""
    if len(symbol) < 2 or not symbol[0].isalpha() or not symbol[1:].isdigit():
        raise ValueError(f"not a valid atom symbol: '{symbol}'")
    form = ord(symbol[0]) - 96
    state = int(symbol[1:])
    if not (0 < form <= 0xFF and 0 <= state <= 0xFF):
        raise ValueError(f"atom symbol out of range: '{symbol}'")
    return (form << 8) | state


def is_bond_species(species: str):
    """'a0.b1' -> ('a0', 'b1'); a plain atom species -> None."""
    if "." in species:
        a, b = species.split(".", 1)
        return a, b
    return None


def _to_int(value, base):
    return int(value, base) if isinstance(value, str) else int(value)


class Reaction:
    def __init__(self, s1, s2, p1, p2, reaction_type, speed, tags=None):
        self.s1, self.s2, self.p1, self.p2 = s1, s2, p1, p2
        self.reaction_type = reaction_type
        self.speed = speed
        self.tags = tags

    @property
    def is_production(self):
        return self.reaction_type == PRODUCTION

    @property
    def is_degradation(self):
        return self.reaction_type == DEGRADATION

    @property
    def is_substrate_contact(self):
        return bool(self.reaction_type & SUBSTRATE_CONTACT)

    @property
    def is_product_contact(self):
        return bool(self.reaction_type & PRODUCT_CONTACT)

    def type_label(self):
        if self.is_production:
            return "production"
        if self.is_degradation:
            return "degradation"
        if self.is_substrate_contact:
            return "contact" if self.is_product_contact else "contact->bond"
        return "link->split" if self.is_product_contact else "link->state"

    def to_string(self):
        """Mirrors Reaction::to_string() in src/reactions.rs."""
        if self.is_production:
            return f"0->{atom_symbol(self.s1)}"
        if self.is_degradation:
            return f"{atom_symbol(self.s1)}->0"
        sub_sep = "+" if self.is_substrate_contact else "."
        prod_sep = "+" if self.is_product_contact else "."
        return (f"{atom_symbol(self.s1)}{sub_sep}{atom_symbol(self.s2)}->"
                f"{atom_symbol(self.p1)}{prod_sep}{atom_symbol(self.p2)}")

    def reactants(self):
        """Species that must already exist for this reaction to be able to fire."""
        if self.is_production:
            return [ENV]
        if self.is_degradation:
            return [atom_symbol(self.s1)]
        if self.is_substrate_contact:
            return [atom_symbol(self.s1), atom_symbol(self.s2)]
        return [bond_symbol(self.s1, self.s2)]  # link reaction: acts on an existing bond

    def products(self):
        """Species this reaction can add, given its reactants are met."""
        if self.is_production:
            return [atom_symbol(self.s1)]
        if self.is_degradation:
            return [ENV]
        if self.is_product_contact:
            return [atom_symbol(self.p1), atom_symbol(self.p2)]  # ends unbonded
        return [bond_symbol(self.p1, self.p2)]  # ends bonded


def load_reactions(path):
    """Load + expand wildcard reactions, mirroring Reaction::news_from_json."""
    with open(path) as f:
        raw = json.load(f)

    reactions = []
    for js in raw:
        s1, s2 = _to_int(js["s1"], 16), _to_int(js["s2"], 16)
        p1, p2 = _to_int(js["p1"], 16), _to_int(js["p2"], 16)
        rtype = _to_int(js["reaction_type"], 2)
        speed = float(js["reaction_speed"])
        tags = js.get("tags")

        s1_wild = (s1 & 0xFF00) == 0xFF00
        s2_wild_all  = (s2 & 0xFF00) == 0xFF00
        s2_wild_same = (s2 & 0xEF00) == 0xEF00

        if s1_wild and s2_wild_all:
            # every (form_a, form_b) combination, states kept from the template
            for fa in range(6):
                for fb in range(16):
                    reactions.append(Reaction(
                        (fa << 8) | (s1 & 0xFF), (fb << 8) | (s2 & 0xFF),
                        (fa << 8) | (p1 & 0xFF), (fb << 8) | (p2 & 0xFF),
                        rtype, speed, tags))
        elif s1_wild and s2_wild_same:
            # same form on both sides, states kept from the template
            for form in range(16):
                reactions.append(Reaction(
                    (form << 8) | (s1 & 0xFF), (form << 8) | (s2 & 0xFF),
                    (form << 8) | (p1 & 0xFF), (form << 8) | (p2 & 0xFF),
                    rtype, speed, tags))
        else:
            reactions.append(Reaction(s1, s2, p1, p2, rtype, speed, tags))

    return reactions


def compute_scope(reactions, initial_species):
    """Fixed-point closure: species reachable from `initial_species` by firing
    any reaction whose reactants are already in scope. Assumes an unbounded
    supply of anything already present — no counts, no kinetics."""
    scope = set(initial_species) | {ENV}
    changed = True
    while changed:
        changed = False
        for r in reactions:
            if all(sp in scope for sp in r.reactants()):
                for sp in r.products():
                    if sp not in scope:
                        scope.add(sp)
                        changed = True
    scope.discard(ENV)
    return scope


def parse_species_token(token):
    """'0x0100' -> 'a0' (hex atom key); anything else passed through as-is
    (symbolic atom 'a0' or bond 'a0.b1')."""
    token = token.strip()
    if token.lower().startswith("0x"):
        return atom_symbol(int(token, 16))
    return token
