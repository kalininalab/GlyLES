from enum import Enum

import networkx as nx

from glyles.glycans.monomer import Monomer


class Chirality(Enum):
    """
    Representation of the chirality of certain atoms based on the Haworth notations of glycans.
    """
    UP = 1
    DOWN = 2
    NONE = 3

    @staticmethod
    def from_string(c):
        if c.lower() == "a":
            return Chirality.DOWN
        if c.lower() == "b":
            return Chirality.UP
        return Chirality.NONE

    def get_opposite(self):
        if self == Chirality.DOWN:
            return Chirality.UP
        if self == Chirality.UP:
            return Chirality.DOWN
        return Chirality.NONE


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


class NXMonomer(Monomer):
    """
    Different monosaccharides that can already be used in the library
    They contain their name, one SMILES representation and a field to store the structure of the glycan once computed
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.__name = kwargs["name"]
        self.__smiles = kwargs["smiles"]
        self.__structure = kwargs["struct"]

    def structure(self):
        """
        Compute and save the structure of this glycan. So far its hard coded on the 5 given types of glycans
        # TODO: Implement a parser for SMILES -> smiles.smiles.SMILES.read()
        Chirality refers to the OH/CH2OH group of the chiral C atom

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
        if self.__structure is None:
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
            if self.__name == "Glc":
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {2: {"chiral": Chirality.DOWN}, 3: {"chiral": Chirality.UP},
                                           4: {"chiral": Chirality.DOWN}, 5: {"chiral": Chirality.UP}})
            elif self.__name == "Fru":
                g.add_edges_from([
                    (10, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 2), (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {3: {"chiral": Chirality.UP}, 4: {"chiral": Chirality.DOWN},
                                           5: {"chiral": Chirality.UP}})
            elif self.__name == "Man":
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15)
                ])
                nx.set_node_attributes(g, {2: {"chiral": Chirality.UP}, 3: {"chiral": Chirality.UP},
                                           4: {"chiral": Chirality.DOWN}, 5: {"chiral": Chirality.UP}})
            elif self.__name == "Gal":
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {2: {"chiral": Chirality.DOWN}, 3: {"chiral": Chirality.UP},
                                           4: {"chiral": Chirality.UP}, 5: {"chiral": Chirality.UP}})
            elif self.__name == "Tal":
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {2: {"chiral": Chirality.UP}, 3: {"chiral": Chirality.UP},
                                           4: {"chiral": Chirality.UP}, 5: {"chiral": Chirality.UP}})
            self.__structure = g

        return self.__structure.copy()

    @staticmethod
    def from_string(mono):
        """
        Get an instance of a glycan according to the provided string representation of the glycan

        Args:
            mono (str): string representation of the glycan of interest

        Returns:
            Glycan according to the monosaccharide provided via mono
        """
        return NXMonomer.__monomers[mono.upper()]


NXMonomer.__monomers = {
    "GLC": NXMonomer(name="Glc", smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", struct=None),
    "FRU": NXMonomer(name="Fru", smiles="C([C@@H]1[C@H]([C@@H](C(O1)(CO)O)O)O)O", struct=None),
    "MAN": NXMonomer(name="Man", smiles="C([C@@H]1[C@H]([C@@H]([C@@H](C(O1)O)O)O)O)O", struct=None),
    "GAL": NXMonomer(name="Gal", smiles="C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O", struct=None),
    "TAL": NXMonomer(name="Tal", smiles="C([C@@H]1[C@@H]([C@@H]([C@@H](C(O1)O)O)O)O)O", struct=None),
}


class DFS:
    """
    Class computing a depth-first search tree from a graph.
    """

    def dfs_tree(self, g, source):
        """
        Compute a depth-first search tree on g starting in the given source node.

        Args:
            g (networkx.Graph): graph to compute the DFS tree for
            source (int): ID of the node to compute the DFS tree from

        Returns:
            A networkx.DiGraph representing the DFS tree of the input graph from the provided source node
        """
        dfst = nx.DiGraph()
        self.__dfs(g, dfst, None, source)

        return dfst

    def __dfs(self, g, dfst, parent, node):
        """
        Recursively perform the depth-first search in the tree. The decision in which order to explore the children is
        important to keep the SMILES-generating algorithm as simple as possible. Here, the search is oriented along the
        ring in the glycan and all side-chains are explored later.

        Args:
            g (networkx.Graph): the graph to compute the DFS tree for
            dfst (networkx.DiGraph): current state of the DFS tree
            parent (int): parent node the search comes from
            node (int): current node id

        Returns:
            Nothing
        """
        # add the current node to the tree
        dfst.add_node(node)

        # in case this node has a parent, add the according edge
        if parent is not None:
            dfst.add_edge(parent, node)

        # iterate over all children in a sorted way, sorted by the membership to the ring and the ID of the node
        children = sorted([c for _, c in list(g.edges(node)) if c != parent],
                          key=lambda c: (not g.nodes[c]["ring"], (c - node) % 10))
        for child in children:
            # if the node is new, explore it
            if child not in dfst.nodes:
                self.__dfs(g, dfst, node, child)

            # if the child already has been seen and its not the parent, we found the ring-closing atoms and mark them
            elif g.nodes[child]["ring"] and g.nodes[node]["ring"]:
                nx.set_node_attributes(dfst, {child: True, node: True}, "ring")


class SMILES:
    """
    Class to hold all methods related to converting a structure into SMILES representation and to generate a structure
    from a SMILES representation.
    """

    def __init__(self):
        self.ring_index = 1

    def write(self, g, source, ring_index=1):
        """
        Compute the SMILES string for the provided graph representing a molecule.
        This only works if called from an Oxygen atom, i.e. a non-chiral atom!

        Args:
            g (networkx.Graph): Graph representation of a molecule
            source (int): Atom-ID to start the SMILES from
            ring_index (int): index of the ring in this specific monosaccharide

        Returns:
            SMILES string of the molecule starting at the given atom
        """
        # compute the DFS tree for the molecule from that given root node
        dfs_tree = DFS().dfs_tree(g, source)

        # based on that DFS tree construct the SMILES
        return self.__construct(g, dfs_tree, source, ring_index)

    def read(self):
        """
        See the exact implementation guidelines in glycans.glycans.Glycan.structure()
        Returns:
            Molecule represented in a graph
        """
        pass

    def __construct(self, g, tree, node, ring_index=1):
        """
        Recursively compute the SMILES from the structure and the according DFS tree. The current step of the search
        has to start in the provided node.
        TODO: Test from non-circular starting points
        Args:
            g (nx.Graph): graph of the molecule
            tree (nx.DiGraph): DFS tree on that molecule
            node (int): node id the construction is currently in

        Returns:
            SMILES representation of the subtree of the current node in tree
        """

        # check for chirality of the molecule. If its not chiral, just add the current atom type to the output
        # IMPORTANT: This method assumes that the chiral C-atoms have an hydrogen atom opposed to the OH group
        if g.nodes[node]["chiral"] != Chirality.NONE:
            """
            parent = next(tree.predecessors(node))
            if parent < 11:
                output = "[C{value}H]".format(value=("@@" if g.nodes[node]["chiral"] == Chirality.DOWN else "@"))
            else:
                output = "[C{value}H]".format(value=("@@" if g.nodes[node]["chiral"] != Chirality.DOWN else "@"))
            """
            # Previously:
            output = "[C{value}H]".format(value=("@@" if g.nodes[node]["chiral"] == Chirality.DOWN else "@"))
        else:
            output = g.nodes[node]["type"].value

        # check if the molecule is start or end of the ring, mark this, remember, there is only one ring in glycans
        if node in nx.get_node_attributes(tree, "ring"):
            output += str(ring_index)

        # recursively step down to the children
        children = list(tree.edges(node))
        if len(children) == 0:  # leaf
            return output
        else:  # children
            for c in [x[1] for x in children]:
                output += "({smiles})".format(smiles=self.__construct(g, tree, c))
            return output
