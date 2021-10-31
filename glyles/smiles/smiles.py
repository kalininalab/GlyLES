class Merger:
    """
    Merge the tree of monomers into a SMILES representation of the complete molecule.
    """

    def __init__(self):
        self.ring_index = 1

    def get_ring_index(self):
        """
        Get index for next ring in the polymer to have disjoint ring indexes in the final SMILES
        
        Returns:
            index for next ring in SMILES
        """
        self.ring_index += 1
        return self.ring_index - 1

    def merge(self, t, root_orientation="n"):
        """
        Merge the provided tree of monomers enriched with the glycans in the nodes and information on the bindings
        between two monomer-nodes in the edges.

        Args:
            t (networkx.DiGraph): Graph representing the glycan to compute the whole SMILES representation for.
            root_orientation (str): Orientation of the root monomer in glycan (can be a (= alpha) or b (= beta))

        Returns:
            SMILES representation as string
        """
        # first mark the atoms that will be replaced in a binding of two monomers
        self.__mark(t, 0, "({}1-?)".format(root_orientation))

        # return the string that can be computed from connecting the monomers as marked above
        return self.__merge(t, 0, 10)

    def __mark(self, t, node, p_edge):
        """
        Recursively mark in every node of the molecule which atoms are being replaced by bound monomers.

        Args:
            t (networkx.DiGraph): Tree representing the glycan to compute the whole SMILES representation for.
            node (int): ID of the node to work on in this method
            p_edge (str): edge annotation to parent monomer

        Returns:
            Nothing
        """
        # get children nodes
        children = [x[1] for x in t.edges(node)]

        # set chirality of atom binding parent
        if p_edge is not None and t.nodes[node]["type"].is_non_chiral():
            print("Chirality-Type:", p_edge[1])
            print(t.nodes[node]["type"].to_smiles())
            t.nodes[node]["type"] = t.nodes[node]["type"].to_chirality(p_edge[1])
            print(t.nodes[node]["type"].to_smiles())
            print(t.nodes[node]["type"].to_chirality(p_edge[1]).to_smiles())

        # check for validity of the tree, i.e. if its a leaf (return, nothing to do) or has too many children (Error)
        if len(children) == 0:  # leaf
            return
        if len(children) > 3:  # too many children
            raise NotImplementedError("Glycans with maximal branching 4 factor not implemented.")

        # iterate over the children and the atoms used to mark binding atoms in my structure
        for child, atom in zip(children, t.nodes[node]["type"].get_dummy_atoms()[0]):
            binding = t.get_edge_data(node, child)["type"]

            t.nodes[node]["type"].mark(int(binding[4]), atom)
            self.__mark(t, child, binding)

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
        me = t.nodes[node]["type"].to_smiles(start, self.get_ring_index())

        # check for validity of the tree, i.e. if its a leaf (return, nothing to do) or has too many children (Error)
        if len(children) == 0:  # leaf
            return me
        if len(children) > 3:  # too many children
            raise NotImplementedError("Glycans with maximal branching factor 4 not implemented.")

        # iterate over the children and the atoms used to mark binding atoms
        for child, atom in zip(children, t.nodes[node]["type"].get_dummy_atoms()[1]):
            binding = list(t.get_edge_data(node, child)["type"])

            child_start = t.nodes[node]["type"].root_atom_id(int(binding[2]))
            # get the SMILES of this child and plug it in in the current own SMILES
            child_smiles = self.__merge(t, child, child_start)
            me = me.replace(atom, child_smiles)

        return me
