import sys
from typing import Union

from networkx.algorithms.isomorphism import DiGraphMatcher
import pydot
from antlr4 import InputStream, CommonTokenStream
from rdkit import Chem

from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.poly.merger import Merger
from glyles.glycans.poly.walker import TreeWalker
from glyles.glycans.utils import ParseError
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.grammar.GlycanParser import GlycanParser


def compare_smiles(c, s):
    Chem.Kekulize(c)
    Chem.Kekulize(s)

    ssmiles = Chem.MolToSmiles(s, kekuleSmiles=True)
    csmiles = Chem.MolToSmiles(c, kekuleSmiles=True)
    return csmiles == ssmiles


def recipe_equality(x, y, full=False):
    """

    Note:
        Neu5Ac and NeuAc and NeuAc are not the same

    Args:
        x:
        y:
        full:

    Returns:

    """
    if not full:
        recipe_x, recipe_y = x.get_recipe(), y.get_recipe()
        return recipe_x[list(zip(*recipe_x))[1].index(GlycanLexer.SAC)] == \
               recipe_y[list(zip(*recipe_y))[1].index(GlycanLexer.SAC)]
    """elif len(recipe_x) == len(recipe_y):
        for x in recipe_x:
            if x not in recipe_y:
                return False
        return True
    return False"""
    return compare_smiles(x.get_structure(), y.get_structure())


class Glycan:
    """
    This class is like an interaction with the Parser for the IUPAC representation of the glycan. The grammar for
    glycans is defined using ANTLR (https://www.antlr.org/). From this ANTLR is able to generate lexer and parser that
    fit the defined grammar. Don't touch those files those are auto generated and therefore mostly uncommented.

    The defined grammar discards the last glycan which is used to define the root of the glycan tree. Therefore, the
    resulting abstract syntax trees (AST)s are not intuitive.
    """

    def __init__(self, iupac, root_orientation="n", start=100, tree_only=False, full=True):
        """
        Initialize the glycan from the IUPAC string.

        Args:
            iupac (str): IUPAC string representation of the glycan to represent
            root_orientation (str): orientation of the root monomer in the glycan (choose from 'a', 'b', 'n')
            start (int): ID of the atom to start with in the root monomer when generating the SMILES
            tree_only (bool): Flag indicating to only parse the tree of glycans and not the modifications
            full (bool): Flag indicating that only fully convertible glycans should be returned, i.e. all modifications
                such as 3-Anhydro-[...] are also present in the SMILES
        """
        self.iupac = iupac
        self.parse_tree = None
        self.glycan_smiles = None
        self.root_orientation = root_orientation
        self.start = start
        self.tree_only = tree_only
        self.factory = MonomerFactory()
        self.full = full
        self.tree_full = True
        self.__parse()

    def count(self, glycan: Union[object, str], exact_nodes=False, exact_edges=False):
        if not isinstance(glycan, Glycan):
            glycan = Glycan(glycan, full=False)

        kwargs = {
            "node_match": lambda x, y: recipe_equality(x["type"], y["type"]),
        }
        if exact_nodes:
            kwargs["node_match"] = lambda x, y: recipe_equality(x["type"], y["type"], True)
        if exact_edges:
            kwargs["edge_match"] = lambda e, f: e["type"] == f["type"]
        matcher = DiGraphMatcher(self.parse_tree, glycan.parse_tree, **kwargs)
        return len(list(matcher.subgraph_isomorphisms_iter()))

    def get_smiles(self):
        """
        Request the SMILES string of the parsed molecule.

        Returns:
            Generated SMILES string
        """
        # return an empty SMILES if the output is required to represent all modifications, but it actually wouldn't
        if not self.tree_only and self.tree_full != self.full:
            return ""

        if self.glycan_smiles is None:
            self.glycan_smiles = Merger(self.factory).merge(self.parse_tree, self.root_orientation, start=self.start)
        self.glycan_smiles = self.glycan_smiles.replace("At", "O-")
        return self.glycan_smiles

    def get_tree(self):
        """
        Request the tree parsed from the IUPAC in this instance.

        Returns:
            The parsed tree with the single monomers in the nodes
        """
        return self.parse_tree

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
            graph.add_edge(pydot.Edge(*edge[::-1], label=self.parse_tree.get_edge_data(*edge)["type"]))
        graph.write(output)

    def __parse(self):
        """
        Adapter on the Lexer and Parser generated by ANTLR based on Grammar.g4.

        Returns:
            Nothing
        """
        # catch the prints of antlr to stderr to check if during parsing an error occurred and the glycan is invalid
        log = []

        class Writer(object):
            @staticmethod
            def write(data):
                log.append(data)

        old_err = sys.stderr
        sys.stderr = Writer()

        # parse the remaining structure description following the grammar, also add the dummy characters
        if not isinstance(self.iupac, str):
            raise ParseError("Only string input can be parsed: " + str(self.iupac))
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
        self.parse_tree, self.tree_full = TreeWalker(self.factory, self.tree_only).parse(tree)

        # if the glycan should be parsed immediately, do so
        if not self.tree_only and self.tree_full == self.full:
            self.glycan_smiles = Merger(self.factory).merge(self.parse_tree, self.root_orientation, start=self.start)
