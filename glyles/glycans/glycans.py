from enum import Enum

import networkx as nx


class Chirality(Enum):
    CC = "@"  # in Haworth notation up
    C = "@@"  # in Haworth notation down
    N = "No"


class Atom(Enum):
    N = "N"
    C = "C"
    O = "O"
    X = "X"
    Y = "Y"
    Z = "Z"


class Glycan(Enum):
    GLC = {"name": "Glc", "smiles": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"}  # from c5 to c1
    FRU = {"name": "Fru", "smiles": "C([C@@H]1[C@H]([C@@H](C(O1)(CO)O)O)O)O"}
    MAN = {"name": "Man", "smiles": "C([C@@H]1[C@H]([C@@H]([C@@H](C(O1)O)O)O)O)O"}
    GAL = {"name": "Gal", "smiles": "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O"}
    TAL = {"name": "Tal", "smiles": "C([C@@H]1[C@@H]([C@@H]([C@@H](C(O1)O)O)O)O)O"}

    def structure(self):
        # TODO: Implement a parser for SMILES -> smiles.smiles.SMILES.read()
        g = nx.Graph()
        g.add_nodes_from([
            (1, {"type": Atom.C, "chiral": Chirality.N}),
            (2, {"type": Atom.C, "chiral": Chirality.N}),
            (3, {"type": Atom.C, "chiral": Chirality.N}),
            (4, {"type": Atom.C, "chiral": Chirality.N}),
            (5, {"type": Atom.C, "chiral": Chirality.N}),
            (6, {"type": Atom.C, "chiral": Chirality.N}),
            (10, {"type": Atom.O, "chiral": Chirality.N}),
            (11, {"type": Atom.O, "chiral": Chirality.N}),
            (12, {"type": Atom.O, "chiral": Chirality.N}),
            (13, {"type": Atom.O, "chiral": Chirality.N}),
            (14, {"type": Atom.O, "chiral": Chirality.N}),
            (15, {"type": Atom.O, "chiral": Chirality.N}),
        ])
        if self == Glycan.GLC:
            g.add_edges_from([
                (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
            ])
            nx.set_node_attributes(g, {2: {"chiral": Chirality.C}, 3: {"chiral": Chirality.CC},
                                       4: {"chiral": Chirality.C}, 5: {"chiral": Chirality.CC}})
        elif self == Glycan.FRU:
            g.add_edges_from([
                (10, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 2), (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
            ])
            nx.set_node_attributes(g, {3: {"chiral": Chirality.CC}, 4: {"chiral": Chirality.C},
                                       5: {"chiral": Chirality.CC}})
        elif self == Glycan.MAN:
            g.add_edges_from([
                (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15)
            ])
            nx.set_node_attributes(g, {2: {"chiral": Chirality.CC}, 3: {"chiral": Chirality.CC},
                                       4: {"chiral": Chirality.C}, 5: {"chiral": Chirality.CC}})
        elif self == Glycan.GAL:
            g.add_edges_from([
                (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
            ])
            nx.set_node_attributes(g, {2: {"chiral": Chirality.C}, 3: {"chiral": Chirality.CC},
                                       4: {"chiral": Chirality.CC}, 5: {"chiral": Chirality.CC}})
        elif self == Glycan.TAL:
            g.add_edges_from([
                (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
            ])
            nx.set_node_attributes(g, {2: {"chiral": Chirality.CC}, 3: {"chiral": Chirality.CC},
                                       4: {"chiral": Chirality.CC}, 5: {"chiral": Chirality.CC}})
        return g


def from_string(mono):
    """

    Args:
        mono (str):

    Returns:
        Glycan according to the monosaccharide provided via mono
    """
    return Glycan[mono.upper()]


print(Glycan.GLC.structure())
