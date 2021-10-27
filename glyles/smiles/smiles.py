import networkx as nx

from glyles.glycans.glycans import Chirality, Atom


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
        children = sorted([c for _, c in list(g.edges(node)) if c != parent], key=lambda c: (not g.nodes[c]["ring"], (c - node) % 10))
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

    def write(self, g, source):
        """
        Compute the SMILES string for the provided graph representing a molecule.
        This only works if called from an Oxygen atom, i.e. a non-chiral atom!

        Args:
            g (networkx.Graph): Graph representation of a molecule
            source (int): Atom-ID to start the SMILES from

        Returns:
            SMILES string of the molecule starting at the given atom
        """
        # compute the DFS tree for the molecule from that given root node
        dfs_tree = DFS().dfs_tree(g, source)

        # based on that DFS tree construct the SMILES
        return self.__construct(g, dfs_tree, source)

    def read(self):
        """
        See the exact implementation guidelines in glycans.glycans.Glycan.structure()
        Returns:
            Molecule represented in a graph
        """
        pass

    def __construct(self, g, tree, node):
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
            output = "[C{value}H]".format(value=("@@" if g.nodes[node]["chiral"] == Chirality.DOWN else "@"))
        else:
            output = g.nodes[node]["type"].value

        # check if the molecule is start or end of the ring, mark this, remember, there is only one ring in glycans
        if node in nx.get_node_attributes(tree, "ring"):
            output += "1"

        # recursively step down to the childen
        children = list(tree.edges(node))
        if len(children) == 0:  # leaf
            return output
        else:  # children
            for c in [x[1] for x in children]:
                output += "({smiles})".format(smiles=self.__construct(g, tree, c))
            return output


class Merger:
    """
    Merge the tree of monomers into a SMILES representation of the complete molecule.
    """

    def merge(self, t):
        """
        Merge the provided tree of monomers enriched with the glycans in the nodes and information on the bindings
        between two monomer-nodes in the edges.

        Args:
            t (networkx.DiGraph): Graph representing the glycan to compute the whole SMILES representation for.

        Returns:
            SMILES representation as string
        """
        # first mark the atoms that will be replaced in a binding of two monomers
        self.__mark(t, 0)

        # return the string that can be computed from connecting the monomers as marked above
        return self.__merge(t, 0, 10)

    def __mark(self, t, node):
        """
        Recursively mark in every node of the molecule which atoms are being replaced by bound monomers.

        Args:
            t (networkx.DiGraph): Graph representing the glycan to compute the whole SMILES representation for.
            node (int): ID of the node to work on in this method

        Returns:
            Nothing
        """
        # get children nodes
        children = [x[1] for x in t.edges(node)]

        # check for validity of the tree, i.e. if its a leaf (return, nothing to do) or has too many children (Error)
        if len(children) == 0:  # leaf
            return
        if len(children) > 3:  # too many children
            raise NotImplementedError("Glycans with maximal branching 4 factor not implemented.")

        # iterate over the children and the atoms used to mark binding atoms
        for child, atom in zip(children, [Atom.X, Atom.Y, Atom.Z]):
            binding = list(t.get_edge_data(node, child)["type"])
            monomer = t.nodes[node]["type"].structure()

            # and find that oxygen atom that will react with an hydrogen when connecting the glycans' monomers
            for _, x in monomer.edges(int(binding[4])):
                if monomer.nodes[x]["type"] == Atom.O and 11 <= x <= 15:
                    monomer.nodes[x]["type"] = atom
                    break

    def __merge(self, t, node, start):
        """
        Recursively merge every node of the molecule with its children and get the SMILES representation of the subtree.

        Args:
            t (networkx.DiGraph): Graph representing the glycan to compute the whole SMILES representation for.
            node (int): ID of the node to work on in this method
            start (int): ID of the atom in the inner graph to start from when generating the SMILES string

        Returns:
            SMILES representation of the subtree of the given node
        """
        # get my children and compute my SMILES string
        children = [x[1] for x in t.edges(node)]
        me = SMILES().write(t.nodes[node]["type"].structure(), start)

        # check for validity of the tree, i.e. if its a leaf (return, nothing to do) or has too many children (Error)
        if len(children) == 0:  # leaf
            return me
        if len(children) > 3:  # too many children
            raise NotImplementedError("Glycans with maximal branching factor 4 not implemented.")

        # iterate over the children and the atoms used to mark binding atoms
        for child, atom in zip(children, ["X", "Y", "Z"]):
            binding = list(t.get_edge_data(node, child)["type"])
            monomer = t.nodes[node]["type"].structure()
            child_start = -1

            # find the ID of the oxygen atom to start the SMILES representation from
            for _, x in monomer.edges(int(binding[2])):
                if monomer.nodes[x]["type"] == Atom.O and 11 <= x <= 15:
                    child_start = x
                    break

            # if the start-point couldn't be found, raise an error
            if child_start == -1:
                raise ValueError("SMILES cannot computed from atom -1!")

            # get the SMILES of this child and plug it in in the current own SMILES
            child_smiles = self.__merge(t, child, child_start)
            me = me.replace(atom, child_smiles)

        return me
