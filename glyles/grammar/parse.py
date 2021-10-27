import os

import networkx as nx
from antlr4 import *
from antlr4.tree.Tree import ParseTree

from glyles.glycans.glycans import Glycan
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.grammar.GlycanParser import GlycanParser

'''
This file is like an interaction with the Parsing of the IUPAC representation of the glycans. The grammar for glycans 
is defined using ANTLR (https://www.antlr.org/). From this ANTLR is able to generate lexer and parser that fit the 
defined grammar. Don't touch those files those are auto generated and therefore mostly uncommented.

The defined grammar discards the last glycan which is used to define the root of the glycan tree. Therefore the 
resulting abstract syntax trees (AST)s are not intuitive.
'''


class TreeWalker:
    def __init__(self):
        """
        Just initialize an empty tree
        """
        self.g = nx.DiGraph()
        self.node_id = 0

    def parse(self, t: ParseTree, init: str):
        """
        Parse an parsed tree (AST) from ANTLR into this networkx graph

        Args:
            t (antlr.ParseTree): result of the parsing step from ANTLR
            init: root monomer of the glycan encoded as string according to the definitions in glycans/glycans.py

        Returns:

        """

        # Initialize the tree with the root
        node_id = self.__add_node(init)

        # add the parsed AST to the networkx graph
        self.__walk(next(t.getChildren()), node_id)

        return self.g

    def __walk(self, t, parent):
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
            node_id = self.__add_node(children[0].symbol.text)
            self.__add_edge(parent, node_id, children[1])
            return node_id

        elif len(children) == 3 and isinstance(children[2], GlycanParser.BranchContext):  # {sac con branch}
            # chain without branching, the parent is the parent of the parsing of the back part
            parent = self.__walk(children[2], parent)
            node_id = self.__add_node(children[0].symbol.text)
            self.__add_edge(parent, node_id, children[1])
            return node_id

        elif len(children) == 3 and isinstance(children[1], GlycanParser.BranchContext):  # {[ branch ]}
            # branching, hand the parent on to the next level
            self.__walk(children[1], parent)
            return parent

        elif len(children) == 7:  # {sac con [ branch ] sac con}
            # branching in a chain, append the end to the parent and hang both branches on that
            node_id = self.__add_node(children[5].symbol.text)
            self.__add_edge(parent, node_id, children[6])
            self.__walk(children[3], node_id)
            node_id2 = self.__add_node(children[0].symbol.text)
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
        self.g.add_edge(parent, child, type=con.symbol.text)

    def __add_node(self, name):
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
        self.g.add_node(node_id, type=Glycan.from_string(name), structure=Glycan.from_string(name).structure())

        return node_id


def parse(smiles, **kwargs):
    """
    Adapter on the Lexer and Parser generated by ANTLR based on Grammar.g4.

    Args:
        smiles (str): IUPAC string of a glycan to be converted to SMILES
        **kwargs: additional arguments

    Returns:
        Tree representation of the glycan with the individual monomers in the nodes of the tree and the binding types
        in the edges
    """

    # Check for monosaccharides and eventually return a single-node-graph
    if "(" not in smiles:
        g = nx.DiGraph()
        g.add_node(0, type=Glycan.from_string(smiles))
        return g

    # Cut off the last monosaccharide as its the root
    bracket_index = max(
        smiles.rindex(")") if ")" in smiles else -1,
        smiles.rindex("]") if "]" in smiles else -1,
    ) + 1
    init = smiles[bracket_index:]
    smiles = smiles[:bracket_index]

    # parse the remaining structure description following the grammar.
    stream = InputStream(data=smiles)
    lexer = GlycanLexer(stream)
    token = CommonTokenStream(lexer)
    parser = GlycanParser(token)
    tree = parser.start()

    # walk through the AST and parse the AST into a networkx representation of the glycan.
    walker = TreeWalker()
    g = walker.parse(tree, init)

    if "output" in kwargs:
        nx.drawing.nx_pydot.write_dot(g, os.path.join(kwargs["output"], f"{smiles}{init}.dot"))

    return g
