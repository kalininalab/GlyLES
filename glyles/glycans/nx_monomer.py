from enum import Enum

import networkx as nx

from glyles.glycans.monomer import Monomer


class NXMonomer(Monomer):
    """
    Different monosaccharides that can already be used in the library.
    They contain their name, one SMILES representation and a field to store the structure of the glycan once computed
    """

    class Chirality(Enum):
        """
        Representation of the chirality of certain atoms based on the Haworth notations of glycans.
        """
        UP = 1
        DOWN = 2
        NONE = 3

        @staticmethod
        def from_string(c):
            """
            Get a conformation flag from a string representation
            Args:
                c (str): name of the conformation in alpha or beta

            Returns:
                Conformation flag according to the char
            """
            if c.lower() == "a":
                return NXMonomer.Chirality.DOWN
            if c.lower() == "b":
                return NXMonomer.Chirality.UP
            return NXMonomer.Chirality.NONE

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

    def __init__(self, origin=None, **kwargs):
        """
        Initialize this class from parent class

        Args:
            origin (Monomer): other instance of Monomers to copy to this instance
            **kwargs: arguments to alternative initialize the monomer if origin is None
        """
        super(NXMonomer, self).__init__(origin, **kwargs)

    def get_dummy_atoms(self):
        """
        Specify some dummy atoms that are used to mark oxygen atoms that will participate in bindings between glycans.
        Here, the atoms will be replaced by instances of the atom enum that are used to define the type of the atom in
        the nodes of the networkx representation of the monomer molecules.

        Returns:
            Two lists:
                * one with atoms that are given as the "atom" argument to the mark-method to replace the oxygen atoms
                * the string representation of the atoms from above, i.e. how the atoms above will be represented in a
                  SMILES string
        """
        return [NXMonomer.Atom.X, NXMonomer.Atom.Y, NXMonomer.Atom.Z], ["X", "Y", "Z"]

    def root_atom_id(self, binding_c_id):
        """
        Get ID of atom that will bind the parent monomer in the glycan. This ID will be given as root argument to
        the to_smiles method.

        Args:
            binding_c_id (int): Integer at which c-position this monomer binds its parent

        Returns:
            id of the atom that binds to the parent, -1 if the root cannot be found
        """
        child_start = -1

        # find the ID of the oxygen atom to start the SMILES representation from
        for _, x in self._get_structure().edges(binding_c_id):
            if self._get_structure().nodes[x]["type"] == NXMonomer.Atom.O and 11 <= x <= 15:
                child_start = x
                break

        # if the start-point couldn't be found, raise an error
        if child_start == -1:
            raise ValueError("SMILES cannot computed from atom -1!")

        return child_start

    def mark(self, position, atom):
        """
        Mark the oxygen atom linked to the carbon atom at the given position ready to participate in the bounding.
        Marking here works based on replacing the oxygen-group bound to the carbon atom at the given position with the
        atom enum instance also provided in the arguments.

        Args:
            position (int): id of the carbon atom whose oxygen atom will from the binding
            atom (object): atom to replace the binding oxygen with

        Returns:
            Nothing
        """
        monomer = self._get_structure()

        # and find that oxygen atom that will react with a hydrogen when connecting the glycans' monomers
        for _, x in monomer.edges(position):
            if monomer.nodes[x]["type"] == NXMonomer.Atom.O and 11 <= x <= 15:
                monomer.nodes[x]["type"] = atom
                break

    def to_smiles(self, root=10, ring_index=1):
        """
        Convert this monomer into a SMILES string representation.
        Use the implementation of the SMILES algorithm fitted to the needs of glycans.

        Args:
            root (int): index of the root atom
            ring_index (int): index of the rings in the atom

        Returns:
            SMILES string representation of this molecule
        """
        return SMILES().write(self._get_structure(), root, ring_index)

    def _get_structure(self):
        """
        Compute and save the structure of this glycan. So far its hard coded on the 4 given types of glycans
        Chirality refers to the OH/CH2OH group of the chiral C atom

        Returns:
            networkx graph representing the structure of the glycan as a graph of its non-hydrogen atoms.
        """
        if self._structure is None:
            g = nx.Graph()
            g.add_nodes_from([
                (1, {"type": NXMonomer.Atom.C, "chiral": NXMonomer.Chirality.NONE, "ring": True}),
                (2, {"type": NXMonomer.Atom.C, "chiral": NXMonomer.Chirality.NONE, "ring": True}),
                (3, {"type": NXMonomer.Atom.C, "chiral": NXMonomer.Chirality.NONE, "ring": True}),
                (4, {"type": NXMonomer.Atom.C, "chiral": NXMonomer.Chirality.NONE, "ring": True}),
                (5, {"type": NXMonomer.Atom.C, "chiral": NXMonomer.Chirality.NONE, "ring": True}),
                (6, {"type": NXMonomer.Atom.C, "chiral": NXMonomer.Chirality.NONE, "ring": False}),
                (10, {"type": NXMonomer.Atom.O, "chiral": NXMonomer.Chirality.NONE, "ring": True}),
                (11, {"type": NXMonomer.Atom.O, "chiral": NXMonomer.Chirality.NONE, "ring": False}),
                (12, {"type": NXMonomer.Atom.O, "chiral": NXMonomer.Chirality.NONE, "ring": False}),
                (13, {"type": NXMonomer.Atom.O, "chiral": NXMonomer.Chirality.NONE, "ring": False}),
                (14, {"type": NXMonomer.Atom.O, "chiral": NXMonomer.Chirality.NONE, "ring": False}),
                (15, {"type": NXMonomer.Atom.O, "chiral": NXMonomer.Chirality.NONE, "ring": False}),
            ])
            if self._name == "Glc":
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {2: {"chiral": NXMonomer.Chirality.DOWN},
                                           3: {"chiral": NXMonomer.Chirality.UP},
                                           4: {"chiral": NXMonomer.Chirality.DOWN},
                                           5: {"chiral": NXMonomer.Chirality.UP}})
            elif self._name == "Man":
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15)
                ])
                nx.set_node_attributes(g, {2: {"chiral": NXMonomer.Chirality.UP},
                                           3: {"chiral": NXMonomer.Chirality.UP},
                                           4: {"chiral": NXMonomer.Chirality.DOWN},
                                           5: {"chiral": NXMonomer.Chirality.UP}})
            elif self._name == "Gal":
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {2: {"chiral": NXMonomer.Chirality.DOWN},
                                           3: {"chiral": NXMonomer.Chirality.UP},
                                           4: {"chiral": NXMonomer.Chirality.UP},
                                           5: {"chiral": NXMonomer.Chirality.UP}})
            elif self._name == "Tal":
                g.add_edges_from([
                    (10, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 10),
                    (1, 11), (2, 12), (3, 13), (4, 14), (5, 6), (6, 15),
                ])
                nx.set_node_attributes(g, {2: {"chiral": NXMonomer.Chirality.UP},
                                           3: {"chiral": NXMonomer.Chirality.UP},
                                           4: {"chiral": NXMonomer.Chirality.UP},
                                           5: {"chiral": NXMonomer.Chirality.UP}})

            if self._config == NXMonomer.Config.ALPHA:
                nx.set_node_attributes(g, {1: {"chiral": NXMonomer.Chirality.DOWN}})
            elif self._config == NXMonomer.Config.BETA:
                nx.set_node_attributes(g, {1: {"chiral": NXMonomer.Chirality.UP}})

            self._structure = g
        return self._structure


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

            # if the child already has been seen, and it's not the parent, we found the ring-closing atoms and mark them
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
        See the exact implementation guidelines in monomer.py description of the monomer description.

        Returns:
            Molecule represented in a graph
        """
        pass

    def __construct(self, g, tree, node, ring_index=1):
        """
        Recursively compute the SMILES from the structure and the according DFS tree. The current step of the search
        has to start in the provided node.

        Args:
            g (nx.Graph): graph of the molecule
            tree (nx.DiGraph): DFS tree on that molecule
            node (int): node id the construction is currently in

        Returns:
            SMILES representation of the subtree of the current node in tree
        """

        # check for chirality of the molecule. If it's not chiral, just add the current atom type to the output
        # IMPORTANT: This method assumes that the chiral C-atoms have a hydrogen atom opposed to the OH group
        if g.nodes[node]["chiral"] != NXMonomer.Chirality.NONE:
            output = "[C{value}H]".format(value=("@@" if g.nodes[node]["chiral"] == NXMonomer.Chirality.DOWN else "@"))
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
