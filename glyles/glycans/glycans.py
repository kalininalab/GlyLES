from enum import Enum

import networkx as nx


class Chirality(Enum):
    """
    Representation of the chirality of certain atoms based on the Haworth notations of glycans.
    """
    UP = 1
    DOWN = 2
    NONE = 3


class Atom(Enum):
    """
    Different types of atoms in the glycans.
    X, Y, Z are used to perform the concatenation of the glycans
    """
    N = "N"
    C = "C"
    O = "O"
    X = "X"
    Y = "Y"
    Z = "Z"


class Glycan(Enum):
    """
    Different monosaccharides that can already be used in the library
    They contain their name, one SMILES representation and a field to store the structure of the glycan once computed
    """
    GLC = {"name": "Glc", "smiles": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", "struct": None}
    FRU = {"name": "Fru", "smiles": "C([C@@H]1[C@H]([C@@H](C(O1)(CO)O)O)O)O", "struct": None}
    MAN = {"name": "Man", "smiles": "C([C@@H]1[C@H]([C@@H]([C@@H](C(O1)O)O)O)O)O", "struct": None}
    GAL = {"name": "Gal", "smiles": "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O", "struct": None}
    TAL = {"name": "Tal", "smiles": "C([C@@H]1[C@@H]([C@@H]([C@@H](C(O1)O)O)O)O)O", "struct": None}

    def structure(self):
        """
        Compute and save the structure of this glycan. So far its hard coded on the 5 given types of glycans
        # TODO: Implement a parser for SMILES -> smiles.smiles.SMILES.read()

        The resulting graph must have some properties (that are especially important for the correct implementation of
        the parser later on:
          - The C atoms have IDs 1 to 6 according to their C1-C6 naming in an glycan.
          - The O atom closing the ring has id 10, so that all other O-atoms have ID 10 higher than the C atom they are
            connected to.
          - The chirality is initial 0 and will only be set in the carbons with respect to the non-OH-group.
          - Additionally, all members of the ring of the molecule have an according flag set to true.

        Returns:
            networkx graph representing the structure of the glycan as a graph of its non-hydrogen atoms.
        """
        if self.value["struct"] is None:
            g = nx.Graph()
            g.add_nodes_from([
                (1, {"type": Atom.C, "chiral": Chirality.NONE, "ring": True}),
                (2, {"type": Atom.C, "chiral": Chirality.NONE, "ring": True}),
                (3, {"type": Atom.C, "chiral": Chirality.NONE, "ring": True}),
                (4, {"type": Atom.C, "chiral": Chirality.NONE, "ring": True}),
                (5, {"type": Atom.C, "chiral": Chirality.NONE, "ring": True}),
                (6, {"type": Atom.C, "chiral": Chirality.NONE, "ring": False}),
                (10, {"type": Atom.O, "chiral": Chirality.NONE, "ring": True}),
                (11, {"type": Atom.O, "chiral": Chirality.NONE, "ring": False}),
                (12, {"type": Atom.O, "chiral": Chirality.NONE, "ring": False}),
                (13, {"type": Atom.O, "chiral": Chirality.NONE, "ring": False}),
                (14, {"type": Atom.O, "chiral": Chirality.NONE, "ring": False}),
                (15, {"type": Atom.O, "chiral": Chirality.NONE, "ring": False}),
            ])
            if self == Glycan.GLC:
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {2: {"chiral": Chirality.DOWN}, 3: {"chiral": Chirality.UP},
                                           4: {"chiral": Chirality.DOWN}, 5: {"chiral": Chirality.UP}})
            elif self == Glycan.FRU:
                g.add_edges_from([
                    (10, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 2), (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {3: {"chiral": Chirality.UP}, 4: {"chiral": Chirality.DOWN},
                                           5: {"chiral": Chirality.UP}})
            elif self == Glycan.MAN:
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15)
                ])
                nx.set_node_attributes(g, {2: {"chiral": Chirality.UP}, 3: {"chiral": Chirality.UP},
                                           4: {"chiral": Chirality.DOWN}, 5: {"chiral": Chirality.UP}})
            elif self == Glycan.GAL:
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {1: {"chiral": Chirality.DOWN}, 2: {"chiral": Chirality.DOWN}, 3: {"chiral": Chirality.UP},
                                           4: {"chiral": Chirality.UP}, 5: {"chiral": Chirality.UP}})
            elif self == Glycan.TAL:
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {2: {"chiral": Chirality.UP}, 3: {"chiral": Chirality.UP},
                                           4: {"chiral": Chirality.UP}, 5: {"chiral": Chirality.UP}})
            self.value["struct"] = g

        return self.value["struct"].copy()

    @staticmethod
    def from_string(mono):
        """
        Get an instance of a glycan according to the provided string representation of the glycan

        Args:
            mono (str): string representation of the glycan of interest

        Returns:
            Glycan according to the monosaccharide provided via mono
        """
        return Glycan[mono.upper()]
