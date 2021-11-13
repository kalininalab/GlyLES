import os
from enum import Enum

import networkx as nx
from antlr4 import *
from rdkit.Chem import MolFromSmiles

from glyles.glycans.nx_monomer import NXMonomer
from glyles.glycans.rdkit_monomer import RDKitMonomer
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.grammar.GlycanParser import GlycanParser
from glyles.utils.merge import Merger


class Glycan:
    """
    This class is like an interaction with the Parser for the IUPAC representation of the glycan. The grammar for
    glycans is defined using ANTLR (https://www.antlr.org/). From this ANTLR is able to generate lexer and parser that
    fit the defined grammar. Don't touch those files those are auto generated and therefore mostly uncommented.

    The defined grammar discards the last glycan which is used to define the root of the glycan tree. Therefore the
    resulting abstract syntax trees (AST)s are not intuitive.
    """

    class TreeWalker:
        def __init__(self):
            """
            Just initialize an empty tree
            """
            self.g = nx.DiGraph()
            self.node_id = 0

        def parse(self, t, mode):
            """
            Parse an parsed tree (AST) from ANTLR into this networkx graph

            Args:
                t (antlr.ParseTree): result of the parsing step from ANTLR
                init (str): root monomer of the glycan encoded as string according to the definitions in glycans/glycans.py
                mode (Glycan.Mode):

            Returns:
                Tree of parsed glycan with monomers in nodes
            """
            children = list(t.getChildren())
            if len(children) == 1:
                self.__add_node(children[0].symbol.text, mode)
                return self.g
            elif len(children) == 2:
                node_id = self.__add_node(children[1].symbol.text, mode)
                self.__walk(children[0], node_id, mode)
                return self.g
            elif len(children) == 3:
                self.__add_node(children[2].symbol.text + children[0].symbol.text, mode)
                return self.g
            elif len(children) == 4:
                node_id = self.__add_node(children[3].symbol.text + children[1].symbol.text, mode)
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

            elif len(children) == 7:  # {sac con [ branch ] sac con}
                # branching in a chain, append the end to the parent and hang both branches on that
                node_id = self.__add_node(children[5].symbol.text, mode)
                self.__add_edge(parent, node_id, children[6])
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
        self.parse()

    def get_smiles(self):
        """
        Request the SMILES string of the parsed molecule.

        Returns:
            Generated SMILES string
        """
        if self.glycan_smiles is None:
            self.glycan_smiles = Merger().merge(self.parse_tree, self.root_orientation, start=self.start)

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
        nx.drawing.nx_pydot.write_dot(self.parse_tree, os.path.join(output, f"{self.iupac}.dot"))

    def parse(self):
        """
        Adapter on the Lexer and Parser generated by ANTLR based on Grammar.g4.

        Returns:
            Nothing
        """
        '''
        # Check for monosaccharides and eventually return a single-node-graph
        if "(" not in self.iupac:
            g = nx.DiGraph()
            if len(self.iupac.split(" ")) == 2:
                g.add_node(0, type=self.monomer_from_string(self.iupac[-1] + self.iupac[:-2], self.mode))
            else:
                g.add_node(0, type=self.monomer_from_string(self.iupac, self.mode))
            self.parse_tree = g
            self.glycan_smiles = g.nodes[0]["type"].get_smiles()
            return

        if self.iupac[-2] == " ":
            root_orientation = self.iupac[-1]
            iupac = self.iupac[:-2]
        else:
            root_orientation = ""
            iupac = self.iupac

        # Cut off the last monosaccharide as its the root
        bracket_index = max(
            iupac.rindex(")") if ")" in iupac else -1,
            iupac.rindex("]") if "]" in iupac else -1,
        ) + 1
        init = root_orientation + iupac[bracket_index:]
        iupac = iupac[:bracket_index]
        '''
        # parse the remaining structure description following the grammar.
        # stream = InputStream(data=self.iupac.replace(" ", "="))
        stream = InputStream(data=self.iupac)
        lexer = GlycanLexer(stream)
        token = CommonTokenStream(lexer)
        parser = GlycanParser(token)
        tree = parser.start()

        # walk through the AST and parse the AST into a networkx representation of the glycan.
        self.parse_tree = Glycan.TreeWalker().parse(tree, self.mode)

        if self.parse_smiles:
            self.glycan_smiles = Merger().merge(self.parse_tree, self.root_orientation, start=self.start)

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
