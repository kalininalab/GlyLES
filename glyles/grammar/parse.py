import os
import sys
from contextlib import contextmanager
from enum import Enum

import networkx as nx
import pydot
from antlr4 import *
from rdkit.Chem import MolFromSmiles

from glyles.glycans.nx_monomer import NXMonomer
from glyles.glycans.rdkit_monomer import RDKitMonomer
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.grammar.GlycanParser import GlycanParser


@contextmanager
def suppress_stdout():
    """
    Source: https://thesmithfam.org/blog/2012/10/25/temporarily-suppress-console-output-in-python/
    Suppress the output of a part of the program
    Use:
    with suppress_stdout():
        // Put code here
    continue with normal code and output
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


class Glycan:
    """
    This class is like an interaction with the Parser for the IUPAC representation of the glycan. The grammar for
    glycans is defined using ANTLR (https://www.antlr.org/). From this ANTLR is able to generate lexer and parser that
    fit the defined grammar. Don't touch those files those are auto generated and therefore mostly uncommented.

    The defined grammar discards the last glycan which is used to define the root of the glycan tree. Therefore, the
    resulting abstract syntax trees (AST)s are not intuitive.
    """

    class __TreeWalker:
        def __init__(self):
            """
            Just initialize an empty tree
            """
            self.g = nx.DiGraph()
            self.node_id = 0

        def parse(self, t, mode):
            """
            Parse a parsed tree (AST) from ANTLR into this networkx graph

            Args:
                t (antlr.ParseTree): result of the parsing step from ANTLR
                mode (Glycan.Mode):

            Returns:
                Tree of parsed glycan with monomers in nodes
            """

            # parse the initial monomer and the orientation of the root monomer
            children = list(t.getChildren())
            if len(children) == 1:  # SAC
                self.__add_node(children[0].symbol.text, mode)
                return self.g
            elif len(children) == 2:  # branch SAC
                node_id = self.__add_node(children[1].symbol.text, mode)
                self.__walk(children[0], node_id, mode)
                return self.g
            elif len(children) == 3:  # SAC ' ' TYPE
                self.__add_node(children[2].symbol.text + "_" + children[0].symbol.text, mode)
                return self.g
            elif len(children) == 4:  # branch SAC ' ' TYPE
                node_id = self.__add_node(children[3].symbol.text + "_" + children[1].symbol.text, mode)
                self.__walk(children[0], node_id, mode)
                return self.g
            else:
                raise RuntimeError("This branch of the if-statement should be unreachable!")

        def __walk(self, t, parent, mode):
            """
            The heart of this class. This method recursively adds nodes to the graph and connects them according to the
            connections in the input IUPAC.

            Args:
                t (antlr.ParseTree): subtree to parse in this recursive step.
                parent (int): ID of the parent node for this (subtree)

            Returns:
                ID of the parent for the remaining recursive procedure after resolving this subtree
            """

            # security check for some nodes that should not occur, but who knows what ANTLR does (or does not) ;-)
            if isinstance(t, ErrorNode):
                raise NotImplementedError("ErrorNodes not implemented yet! (Should also never occur!)")
            elif isinstance(t, TerminalNode):
                raise NotImplementedError("Terminal nodes should be unreachable!")

            children = list(t.getChildren())
            if len(children) == 2:  # {sac con}
                # terminal element, add the node with the connection
                node_id = self.__add_node(children[0].symbol.text, mode)
                self.__add_edge(parent, node_id, children[1])
                return node_id

            elif len(children) == 3 and isinstance(children[2], GlycanParser.BranchContext):  # {sac con branch}
                # chain without branching, the parent is the parent of the parsing of the back part
                parent = self.__walk(children[2], parent, mode)
                node_id = self.__add_node(children[0].symbol.text, mode)
                self.__add_edge(parent, node_id, children[1])
                return node_id

            elif len(children) == 3 and isinstance(children[1], GlycanParser.BranchContext):  # {[ branch ]}
                # branching, hand the parent on to the next level
                self.__walk(children[1], parent, mode)
                return parent

            elif len(children) == 6:  # {sac con [ branch ] branch}
                # branching in a chain, append the end to the parent and hang both branches on that
                node_id = self.__walk(children[5], parent, mode)
                self.__walk(children[3], node_id, mode)
                node_id2 = self.__add_node(children[0].symbol.text, mode)
                self.__add_edge(node_id, node_id2, children[1])
                return node_id2

            # there should be no case missing, but who knows...
            raise NotImplementedError("This should be unreachable")

        def __add_edge(self, parent, child, con):
            """
            Add an edge between the provided ids of the parent and the children in the glycan tree.

            Args:
                parent (int): ID of parent in the connection of the glycan tree
                child (int): ID of the child in the connection of the glycan tree
                con: Connection element from the ParseTree

            Returns:
                Nothing
            """
            self.g.add_edge(parent, child, type=str(con))

        def __add_node(self, name, mode):
            """
            Add a new node to the network based on the name of the represented glycan.

            Args:
                name (str): Name of the glycan to be represented in the new node

            Returns:
                ID of the newly added node
            """
            # get the id for the node and increase the id for the next node
            node_id = self.node_id
            self.node_id += 1

            # add the node to the network and store the enum of the glycan as attribute
            self.g.add_node(
                node_id,
                type=Glycan.monomer_from_string(name, mode),
            )

            return node_id

    class __Merger:
        """
        Merge the tree of monomers into a SMILES representation of the complete molecule.
        """

        def merge(self, t, root_orientation="n", start=10):
            """
            Merge the provided tree of monomers enriched with the glycans in the nodes and information on the bindings
            between two monomer-nodes in the edges.

            Args:
                t (networkx.DiGraph): Graph representing the glycan to compute the whole SMILES representation for.
                root_orientation (str): Orientation of the root monomer in glycan (can be a (= alpha) or b (= beta))
                start (int): index of the atom to start the SMILES generation from

            Returns:
                SMILES representation as string
            """
            # first mark the atoms that will be replaced in a binding of two monomers
            self.__mark(t, 0, "({}1-?)".format(root_orientation))

            # return the string that can be computed from connecting the monomers as marked above
            return self.__merge(t, 0, start, 1)

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
                t.nodes[node]["type"] = t.nodes[node]["type"].to_chirality(p_edge[1])

            # check for validity of the tree, ie if it's a leaf (return, nothing to do) or has too many children (Error)
            if len(children) == 0:  # leaf
                return
            if len(children) > 3:  # too many children
                raise NotImplementedError("Glycans with maximal branching 4 factor not implemented.")

            # iterate over the children and the atoms used to mark binding atoms in my structure
            for child, atom in zip(children, t.nodes[node]["type"].get_dummy_atoms()[0]):
                binding = t.get_edge_data(node, child)["type"]

                t.nodes[node]["type"].mark(int(binding[4]), atom)
                self.__mark(t, child, binding)

        def __merge(self, t, node, start, ring_index):
            """
            Recursively merge every node of the molecule with its children and get the SMILES representation of the
            subtree.

            Args:
                t (networkx.DiGraph): Graph representing the glycan to compute the whole SMILES representation for.
                node (int): ID of the node to work on in this method
                start (int): ID of the atom in the inner graph to start from when generating the SMILES string

            Returns:
                SMILES representation of the subtree of the given node
            """
            # get my children and compute my SMILES string
            children = [x[1] for x in t.edges(node)]
            me = t.nodes[node]["type"].to_smiles(start, ring_index)

            # check for validity of the tree, ie if it's a leaf
            if len(children) == 0:  # leaf
                return me

            # iterate over the children and the atoms used to mark binding atoms
            for child, atom in zip(children, t.nodes[node]["type"].get_dummy_atoms()[1]):
                binding = list(t.get_edge_data(node, child)["type"])

                child_start = t.nodes[node]["type"].root_atom_id(int(binding[2]))
                if child_start == -1:
                    raise ValueError("No child start found.")

                # get the SMILES of this child and plug it in the current own SMILES
                child_smiles = self.__merge(t, child, child_start, ring_index + 1)
                me = me.replace(atom, child_smiles)
            return me

    class Mode(Enum):
        """
        Enumerate different modes how to represent the monomers in the tree.
        """
        DEFAULT_MODE = "rdkit"
        NETWORKX_MODE = "nx"
        RDKIT_MODE = "rdkit"

    def __init__(self, iupac, mode=Mode.DEFAULT_MODE, root_orientation="n", start=10, parse=True):
        """
        Initialize the glycan from the IUPAC string.

        Args:
            iupac (str): IUPAC string representation of the glycan to represent
            mode (Glycan.Mode): Glycan mode to use. This controls the representation used (either networkx or RDKit).
        """
        self.iupac = iupac
        self.mode = mode
        self.parse_tree = None
        self.glycan_smiles = None
        self.root_orientation = root_orientation
        self.start = start
        self.parse_smiles = parse
        if parse:
            self.__parse()

    def get_smiles(self):
        """
        Request the SMILES string of the parsed molecule.

        Returns:
            Generated SMILES string
        """
        if self.glycan_smiles is None:
            self.glycan_smiles = Glycan.__Merger().merge(self.parse_tree, self.root_orientation, start=self.start)

        return self.glycan_smiles

    def get_tree(self):
        """
        Request the tree parsed from the IUPAC in this instance.

        Returns:
            The parsed tree with the single monomers in the nodes
        """
        return self.parse_tree

    def to_pdb(self, output):
        """
        Compute a 3-dimensional conformation of the molecule using ForceFields. The result will be stored in a PDB file.

        Args:
            output (str): Path to store the PDB file in

        Returns:
            Nothing
        """
        mol = MolFromSmiles(self.glycan_smiles)
        raise NotImplementedError("PDB conversion planned, but not implemented yet.")

    def save_dot(self, output):
        """
        Save the tree structure of the encoded glycan molecule into a dot file visualizing the graph of monomers.

        Args:
            output (str): path to store the DOT file in

        Returns:
            Nothing
        """
        graph = pydot.Dot("iupac_tree")
        for node in range(len(self.parse_tree.nodes)):
            graph.add_node(pydot.Node(node, label=self.parse_tree.nodes[node]["type"].get_name()))
        for edge in self.parse_tree.edges():
            graph.add_edge(pydot.Edge(*edge, label=self.parse_tree.get_edge_data(*edge)["type"]))
        graph.write(output + ".dot")

    def __parse(self):
        """
        Adapter on the Lexer and Parser generated by ANTLR based on Grammar.g4.

        Returns:
            Nothing
        """
        # parse the remaining structure description following the grammar.
        # with suppress_stdout():
        stream = InputStream(data=self.iupac)
        lexer = GlycanLexer(stream)
        token = CommonTokenStream(lexer)
        parser = GlycanParser(token)
        tree = parser.start()

        # walk through the AST and parse the AST into a networkx representation of the glycan.
        self.parse_tree = Glycan.__TreeWalker().parse(tree, self.mode)

        if self.parse_smiles:
            self.glycan_smiles = Glycan.__Merger().merge(self.parse_tree, self.root_orientation, start=self.start)

    @staticmethod
    def monomer_from_string(mono, mode):
        """
        Get a monomer representation from its short name.

        Args:
            mono (str): monomer in three (or more) letter representation, i.e. Glc = Glucose, Gal = Galactose, etc.
            mode (Glycan.Mode): mode determining which monomer-mode to use

        Returns:
            An instance of a monomer representation for this graph
        """
        if mode == Glycan.Mode.NETWORKX_MODE:
            return NXMonomer.from_string(mono)
        elif mode == Glycan.Mode.RDKIT_MODE:
            return RDKitMonomer.from_string(mono)
        else:
            raise ValueError("Unknown representation mode for monomers!")


if __name__ == '__main__':  # testing proposes
    i = 1
    if i == 0:
        print(RDKitMonomer.from_string("glc").to_smiles(-1, 0))
    elif i == 1:
        glycan = Glycan("Man(a1-2)Man", mode=Glycan.Mode.RDKIT_MODE)
        print(glycan.get_smiles())
    print("Done")
