from enum import Enum

import networkx as nx

"""
SMILES:
 @@ : clockwise rotation (in Haworth down)
 @  : counter-clockwise rotation (in Haworth up)
"""


class Atom(Enum):
    N = "N"
    C = "C"
    O = "O"
    X = "X"


class Glycan(Enum):
    GLC = {"name": "Glc", "smiles": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"}  # from c5 to c1
    FRU = {"name": "Fru", "smiles": "C([C@@H]1[C@H]([C@@H](C(O1)(CO)O)O)O)O"}
    MAN = {"name": "Man", "smiles": "C([C@@H]1[C@H]([C@@H]([C@@H](C(O1)O)O)O)O)O"}
    GAL = {"name": "Gal", "smiles": "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O"}
    TAL = {"name": "Tal", "smiles": "C([C@@H]1[C@@H]([C@@H]([C@@H](C(O1)O)O)O)O)O"}

    def structure(self):
        g = nx.Graph()
        g.add_nodes_from([
            (1, {"type": Atom.C}),
            (2, {"type": Atom.C}),
            (3, {"type": Atom.C}),
            (4, {"type": Atom.C}),
            (5, {"type": Atom.C}),
            (6, {"type": Atom.C}),
            (10, {"type": Atom.O}),
            (11, {"type": Atom.O}),
            (12, {"type": Atom.O}),
            (13, {"type": Atom.O}),
            (14, {"type": Atom.O}),
            (15, {"type": Atom.O}),
        ])
        if self == Glycan.GLC:
            g.add_edges_from([
                (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12, {"c": "down"}), (3, 13, {"c": "up"}), (4, 14, {"c": "down"}), (6, 15, {"c": "up"}),
            ])
        elif self == Glycan.FRU:
            g.add_edges_from([
                (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12), (3, 13, {"c": "up"}), (4, 14, {"c": "down"}), (5, 6, {"c": "up"}), (6, 15),
            ])
        elif self == Glycan.MAN:
            g.add_edges_from([
                (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12, {"c": "up"}), (3, 13, {"c": "up"}), (4, 14, {"c": "down"}), (6, 15, {"c": "up"}),
            ])
        elif self == Glycan.GAL:
            g.add_edges_from([
                (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12, {"c": "down"}), (3, 13, {"c": "up"}), (4, 14, {"c": "up"}), (6, 15, {"c": "up"}),
            ])
        elif self == Glycan.TAL:
            g.add_edges_from([
                (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                (1, 11), (2, 12, {"c": "up"}), (3, 13, {"c": "up"}), (4, 14, {"c": "up"}), (6, 15, {"c": "up"}),
            ])
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
