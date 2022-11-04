import copy
import re

import numpy as np

from glyles.glycans.utils import sanitize_smiles


class Merger:
    """
    Merge the tree of monomers into a SMILES representation of the complete molecule.
    """

    def __init__(self, factory):
        """
        Create a merger class to merge a parsed glycan tree into a SMILES string.

        Args:
            factory (MonomerFactory): factory instance to use to generate the monomers for the glycan tree from
        """
        self.factory = factory

    def merge(self, t, root_orientation="n", start=100):
        """
        Merge the provided tree of monomers enriched with the glycans in the nodes and information on the bindings
        between two monomer-nodes in the edges. The input graph is not changed during this process, as the methods
        deep-copies the glycan representing datastructure.

        Args:
            t (networkx.DiGraph): Graph representing the glycan to compute the whole SMILES representation for.
            root_orientation (str): Orientation of the root monomer in glycan (can be a (= alpha) or b (= beta))
            start (int): index of the atom to start the SMILES generation from

        Returns:
            SMILES representation as string
        """
        # deep-copy the tree to not modify the provided argument
        t = copy.deepcopy(t)

        # first mark the atoms that will be replaced in a binding of two monomers
        self.__mark(t, 0, f"({root_orientation}1-?)")

        # return the string that can be computed from connecting the monomers as marked above
        if np.where(t.nodes[0]["type"].get_features()[:, 1] == start)[0].size != 0:
            position = int(np.argwhere(t.nodes[0]["type"].get_features()[:, 1] == start).squeeze())
        else:
            position = int(np.argwhere(t.nodes[0]["type"].get_features()[:, 1] == 1).squeeze())
        return self.__merge(t, 0, position, 0)

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
            t.nodes[node]["type"] = t.nodes[node]["type"].to_chirality(p_edge[1], self.factory)

        # check for validity of the tree, ie if it's a leaf (return, nothing to do) or has too many children (Error)
        if len(children) == 0:  # leaf
            return
        if len(children) > 4:  # too many children
            raise NotImplementedError("Glycans with maximal branching factor greater then 3 not implemented.")

        # iterate over the children and the atoms used to mark binding atoms in my structure
        for child, atom in zip(children, t.nodes[node]["type"].get_dummy_atoms()):
            binding = re.findall(r'\d+', t.get_edge_data(node, child)["type"])[1]

            t.nodes[node]["type"].mark(int(binding), *atom)
            self.__mark(t, child, t.get_edge_data(node, child)["type"])

    def __merge(self, t, node, start, ring_index):
        """
        Recursively merge every node of the molecule with its children and get the SMILES representation of the
        subtree.

        Args:
            t (networkx.DiGraph): Graph representing the glycan to compute the whole SMILES representation for.
            node (int): ID of the node to work on in this method
            start (int): ID of the atom in the inner graph to start from when generating the SMILES string
            ring_index (int): Index of the ring to use when generating the SMILES strings

        Returns:
            SMILES representation of the subtree of the given node
        """
        # get my children and compute my SMILES string
        children = [x[1] for x in t.edges(node)]
        me = t.nodes[node]["type"].to_smiles(ring_index, root_id=start)

        # check for validity of the tree, ie if it's a leaf
        if len(children) == 0:  # leaf
            return me

        # iterate over the children and the atoms used to mark binding atoms
        for child, (o_atom, n_atom) in zip(children, t.nodes[node]["type"].get_dummy_atoms()):
            binding = re.findall(r'\d+', t.get_edge_data(node, child)["type"])[0]

            child_start = t.nodes[child]["type"].root_atom_id(int(binding))
            if child_start == -1:
                raise ValueError("No child start found.")

            # get the SMILES of this child and plug it in the current own SMILES
            child_smiles = self.__merge(t, child, child_start, ring_index + 1)
            if o_atom[1] in me:
                me = re.sub(o_atom[2], child_smiles, me)
            elif n_atom[1] in me:
                me = re.sub(n_atom[2], "N(" + child_smiles[1:] + ")", me)
            me = sanitize_smiles(me)
            # me = me.replace("((", "(").replace("))", ")")
            # me = me.replace(atom, child_smiles)
        return me
