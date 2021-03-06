import sys

import networkx as nx
import numpy as np
import pydot
from antlr4 import *
from rdkit.Chem import MolFromSmiles

from glyles.glycans.utils import Mode, UnreachableError, ParseError
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.grammar.GlycanParser import GlycanParser


class Glycan:
    """
    This class is like an interaction with the Parser for the IUPAC representation of the glycan. The grammar for
    glycans is defined using ANTLR (https://www.antlr.org/). From this ANTLR is able to generate lexer and parser that
    fit the defined grammar. Don't touch those files those are auto generated and therefore mostly uncommented.

    The defined grammar discards the last glycan which is used to define the root of the glycan tree. Therefore, the
    resulting abstract syntax trees (AST)s are not intuitive.
    """

    class __TreeWalker:
        def __init__(self, factory, tree_only):
            """
            Just initialize an empty tree

            Args:
                factory (MonomerFactory): factory instance to use to generate the monomers for the glycan tree from
            tree_only (bool): Flag indicating to only parse the tree of glycans and not the modifications
            """
            self.g = nx.DiGraph()
            self.factory = factory
            self.node_id = 0
            self.tree_only = tree_only

        def parse(self, t, mode):
            """
            Parse a parsed tree (AST) from ANTLR into this networkx graph

            Args:
                t (antlr.ParseTree): result of the parsing step from ANTLR
                mode (Mode): mode determining which monomer-mode to use

            Returns:
                Tree of parsed glycan with monomers in nodes
            """

            # parse the initial monomer and the orientation of the root monomer
            # and remove the first and last char, i.e. '{' and '}'
            children = list(t.getChildren())[1:-1]
            if len(children) == 1:  # glycan
                self.__add_node(children[0], mode)
                return self.g
            elif len(children) == 2:  # branch glycan
                node_id = self.__add_node(children[1], mode)
                self.__walk(children[0], node_id, mode)
                return self.g
            elif len(children) == 3:  # SAC ' ' TYPE
                self.__add_node(children[0], mode, children[2].symbol.text)
                return self.g
            elif len(children) == 4:  # branch SAC ' ' TYPE
                node_id = self.__add_node(children[1], mode, children[3].symbol.text)
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
                mode (Mode): mode determining which monomer-mode to use

            Returns:
                ID of the parent for the remaining recursive procedure after resolving this subtree
            """

            # security check for some nodes that should not occur, but who knows what ANTLR does (or does not) ;-)
            if isinstance(t, ErrorNode):
                raise UnreachableError("ErrorNodes in parsing!")
            elif isinstance(t, TerminalNode):
                raise UnreachableError("TerminalNodes in parsing!")

            children = list(t.getChildren())
            if len(children) == 2:  # {glycan con}
                # terminal element, add the node with the connection
                node_id = self.__add_node(children[0], mode)
                self.__add_edge(parent, node_id, children[1])
                return node_id

            elif len(children) == 3 and isinstance(children[2], GlycanParser.BranchContext):  # {glycan con branch}
                # chain without branching, the parent is the parent of the parsing of the back part
                parent = self.__walk(children[2], parent, mode)
                node_id = self.__add_node(children[0], mode)
                self.__add_edge(parent, node_id, children[1])
                return node_id

            elif len(children) == 3 and isinstance(children[1], GlycanParser.BranchContext):  # {'[' branch ']'}
                # branching, hand the parent on to the next level
                self.__walk(children[1], parent, mode)
                return parent

            elif len(children) == 6:  # {glycan con '[' branch ']' branch}
                # branching in a chain, append the end to the parent and hang both branches on that
                node_id = self.__walk(children[5], parent, mode)
                self.__walk(children[3], node_id, mode)
                node_id2 = self.__add_node(children[0], mode)
                self.__add_edge(node_id, node_id2, children[1])
                return node_id2

            # there should be no case missing, but who knows...
            raise UnreachableError("Invalid branching in glycan tree")

        def __add_edge(self, parent, child, con):
            """
            Add an edge between the provided ids of the parent and the children in the glycan tree.

            Args:
                parent (int): ID of parent in the connection of the glycan tree
                child (int): ID of the child in the connection of the glycan tree
                con (str): Connection element from the ParseTree

            Returns:
                Nothing
            """
            self.g.add_edge(parent, child, type=str(con))

        def __add_node(self, node, mode, config=""):
            """
            Add a new node to the network based on the name of the represented glycan.

            Args:
                node: Node of the parsed tree to be parsed into a monomer
                mode (Mode): mode determining which monomer-mode to use
                config (str): configuration to be applied to the monomer

            Returns:
                ID of the newly added node
            """
            # get the id for the node and increase the id for the next node
            node_id = self.node_id
            self.node_id += 1

            # add the node to the network and store the enum of the glycan as attribute
            recipe = []
            for child in node.children:
                if isinstance(child, GlycanParser.DerivContext):
                    for c in child.children:
                        recipe.append((str(c), c.symbol.type))
                else:
                    recipe.append((str(child), child.symbol.type))
            self.g.add_node(
                node_id,
                type=self.factory.create(recipe, mode, config, tree_only=self.tree_only),
            )

            return node_id

    class __Merger:
        """
        Merge the tree of monomers into a SMILES representation of the complete molecule.
        """

        def __init__(self, factory):
            """
            Create a merger class to merge a parsed glycan tree into a SMILES string

            Args:
                factory (MonomerFactory): factory instance to use to generate the monomers for the glycan tree from
            """
            self.factory = factory

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
            position = int(np.argwhere(t.nodes[0]["type"].get_features()[:, 1] == start).squeeze())
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
            if len(children) > 2:  # too many children
                raise NotImplementedError("Glycans with maximal branching factor greater then 2 not implemented.")

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
                ring_index (int): Index of the ring to use when generating the SMILES strings

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

                # child_start = t.nodes[node]["type"].root_atom_id(int(binding[2]))
                child_start = t.nodes[child]["type"].root_atom_id(int(binding[2]))
                if child_start == -1:
                    raise ValueError("No child start found.")

                # get the SMILES of this child and plug it in the current own SMILES
                child_smiles = self.__merge(t, child, child_start, ring_index + 1)
                me = me.replace(atom, child_smiles)
            return me

    def __init__(self, iupac, factory, mode=Mode.DEFAULT_MODE, root_orientation="n", start=10, tree_only=False):
        """
        Initialize the glycan from the IUPAC string.

        Args:
            iupac (str): IUPAC string representation of the glycan to represent
            factory (MonomerFactory): factory instance to use to generate the monomers for the glycan tree from
            mode (Mode): Glycan mode to use. This controls the representation used (either networkx or RDKit)
            root_orientation (str): orientation of the root monomer in the glycan (choose from 'a', 'b', 'n')
            start (int): ID of the atom to start with in the root monomer when generating the SMILES
            tree_only (bool): Flag indicating to only parse the tree of glycans and not the modifications
        """
        self.iupac = iupac
        self.mode = mode
        self.parse_tree = None
        self.glycan_smiles = None
        self.root_orientation = root_orientation
        self.start = start
        self.tree_only = tree_only
        self.factory = factory
        self.__parse()

    def get_smiles(self):
        """
        Request the SMILES string of the parsed molecule.

        Returns:
            Generated SMILES string
        """
        if self.glycan_smiles is None:
            self.glycan_smiles = Glycan.__Merger(self.factory).merge(self.parse_tree, self.root_orientation,
                                                                     start=self.start)
        self.glycan_smiles = self.glycan_smiles.replace("At", "O-")
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

    def save_dot(self, output, horizontal=False):
        """
        Save the tree structure of the encoded glycan molecule into a dot file visualizing the graph of monomers.

        Args:
            output (str): path to store the DOT file in
            horizontal (bool): Show graph in horizontal orientation from left to right

        Returns:
            Nothing
        """
        if horizontal:
            graph = pydot.Dot("iupac_tree", rankdir="LR")
        else:
            graph = pydot.Dot("iupac_tree")
        for node in range(len(self.parse_tree.nodes)):
            graph.add_node(pydot.Node(node, label=self.parse_tree.nodes[node]["type"].get_name(full=True)))
        for edge in self.parse_tree.edges():
            graph.add_edge(pydot.Edge(*edge, label=self.parse_tree.get_edge_data(*edge)["type"]))
        graph.write(output + ".dot")

    def __parse(self):
        """
        Adapter on the Lexer and Parser generated by ANTLR based on Grammar.g4.

        Returns:
            Nothing
        """
        # catch the prints of antlr to stderr to check if during parsing an error occurred and the glycan is invalid
        log = []

        class writer(object):
            @staticmethod
            def write(data):
                log.append(data)

        old_err = sys.stderr
        sys.stderr = writer()

        # parse the remaining structure description following the grammar, also add the dummy characters
        stream = InputStream(data='{' + self.iupac + '}')
        lexer = GlycanLexer(stream)
        token = CommonTokenStream(lexer)
        parser = GlycanParser(token)
        tree = parser.start()

        sys.stderr = old_err

        # if the glycan is invalid, set its structure to None and the SMILES string to empty and return
        if len(log) != 0:
            self.parse_tree = None
            self.glycan_smiles = ""
            raise ParseError("Glycan cannot be parsed:\n" + log[0])

        # walk through the AST and parse the AST into a networkx representation of the glycan.
        self.parse_tree = Glycan.__TreeWalker(self.factory, self.tree_only).parse(tree, self.mode)

        # if the glycan should be parsed immediately, do so
        if not self.tree_only:
            self.glycan_smiles = Glycan.__Merger(self.factory).merge(self.parse_tree, self.root_orientation,
                                                                     start=self.start)
