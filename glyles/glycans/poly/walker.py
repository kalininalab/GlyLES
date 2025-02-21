import networkx as nx
from antlr4 import ErrorNode, TerminalNode

from glyles.glycans.utils import UnreachableError, ketoses2
from glyles.gwb.GWBParser import GWBParser
from glyles.iupac.IUPACLexer import IUPACLexer
from glyles.iupac.IUPACParser import IUPACParser
from glyles.gwb.GWBLexer import GWBLexer


class TreeWalker:
    def __init__(self, factory, tree_only):
        """
        Just initialize an empty tree.

        Args:
            factory (MonomerFactory): factory instance to use to generate the monomers for the glycan tree from
        tree_only (bool): Flag indicating to only parse the tree of glycans and not the modifications
        """
        self.g = nx.DiGraph()
        self.factory = factory
        self.node_id = 0
        self.tree_only = tree_only
        self.full = True

    def parse(self, t):
        for c in filter(lambda x: not isinstance(x, (TerminalNode, ErrorNode)), t.getChildren()):
            if isinstance(c, IUPACParser.BranchContext):
                self.walk(c, self.node_id)
            elif isinstance(c, IUPACParser.BeginContext):
                self.parse_int(c)
        return self.g, self.full and len(list(nx.connected_components(self.g.to_undirected()))) == 1

    def parse_int(self, t):
        """
        Parse a parsed tree (AST) from ANTLR into this networkx graph.

        Args:
            t (IUPACParser.BeginContext): result of the parsing step from ANTLR

        Returns:
            Tree of parsed glycan with monomers in nodes
        """

        # parse the initial monomer and the orientation of the root monomer
        # and remove the first and last char, i.e. '{' and '}'
        # if isinstance(t, )
        children = list(t.getChildren())
        if len(children) == 1:  # glycan
            self.add_node(children[0])
        elif len(children) == 2:  # branch glycan
            node_id = self.add_node(children[1])
            self.walk(children[0], node_id)
        elif len(children) == 3:  # SAC ' ' TYPE
            self.add_node(children[0], config=children[2].symbol.text)
        elif len(children) == 4:  # branch SAC ' ' TYPE
            node_id = self.add_node(children[1], config=children[3].symbol.text)
            self.walk(children[0], node_id)
        else:
            raise RuntimeError("This branch of the if-statement should be unreachable!")

    def walk(self, t, parent):
        """
        The heart of this class. This method recursively adds nodes to the graph and connects them according to the
        connections in the input IUPAC.

        Args:
            t (IUPACParser.BranchContext): subtree to parse in this recursive step.
            parent (int): ID of the parent node for this (subtree)

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
            node_id = self.add_node(children[0])
            self.full &= self.add_edge(parent, node_id, children[1])
            return node_id

        elif len(children) == 3 and isinstance(children[2], IUPACParser.BranchContext):  # {glycan con branch}
            # chain without branching, the parent is the parent of the parsing of the back part
            parent = self.walk(children[2], parent)
            node_id = self.add_node(children[0])
            self.full &= self.add_edge(parent, node_id, children[1])
            return node_id

        elif len(children) == 3 and isinstance(children[1], IUPACParser.BranchContext):  # {'[' branch ']'}
            # branching, hand the parent on to the next level
            self.walk(children[1], parent)
            return parent

        elif len(children) == 6:  # {glycan con '[' branch ']' branch}
            # branching in a chain, append the end to the parent and hang both branches on that
            node_id = self.walk(children[5], parent)
            self.walk(children[3], node_id)
            node_id2 = self.add_node(children[0])
            self.full &= self.add_edge(node_id, node_id2, children[1])
            return node_id2

        elif len(children) == 9:  # {glycan con '[' branch ']' '[' branch ']' branch}
            # branching in a chain, append the end to the parent and hang both branches on that
            node_id = self.walk(children[8], parent)
            self.walk(children[3], node_id)
            self.walk(children[6], node_id)
            node_id2 = self.add_node(children[0])
            self.full &= self.add_edge(node_id, node_id2, children[1])
            return node_id2

        elif len(children) == 12:  # {glycan con '[' branch ']' '[' branch ']' '[' branch ']' branch}
            # branching in a chain, append the end to the parent and hang both branches on that
            node_id = self.walk(children[11], parent)
            self.walk(children[3], node_id)
            self.walk(children[6], node_id)
            self.walk(children[9], node_id)
            node_id2 = self.add_node(children[0])
            self.full &= self.add_edge(node_id, node_id2, children[1])
            return node_id2

        # there should be no case missing, but who knows...
        raise UnreachableError("Invalid branching in glycan tree")

    def context2str(self, node):
        output = ""
        for c in node.getChildren():
            if isinstance(c, TerminalNode):
                output += str(c)
            else:
                output += self.context2str(c)
        return output

    def add_edge(self, parent, child, con) -> bool:
        """
        Add an edge between the provided ids of the parent and the children in the glycan tree.

        Args:
            parent (int): ID of parent in the connection of the glycan tree
            child (int): ID of the child in the connection of the glycan tree
            con (str): Connection element from the ParseTree

        Returns:
            Nothing
        """
        # for floating elements as they add a new connected component
        if parent == child:
            return True

        con = self.context2str(con)
        if "(" not in con and ")" not in con:
            if "-" not in con:
                bond = ("2-" if (self.g.nodes[child]["type"].get_lactole,
                                 self.g.nodes[child]['type'].get_name()) in ketoses2 else "1-")
                con = con[0] + bond + con[1:]
            con = "(" + con + ")"
        self.g.add_edge(parent, child, type=con)
        # Identify standard and narrow wildcard connections
        return "?" not in con and "/" not in con

    def build_recipe(self, node):
        recipe = []
        for c in node.getChildren():
            if isinstance(c, TerminalNode):
                recipe.append((str(c), IUPACLexer.SAC if isinstance(node, (IUPACParser.SaciContext, GWBParser.SaciContext)) else c.symbol.type))
            else:
                tmp = self.build_recipe(c)
                if isinstance(c, (IUPACParser.ModiContext, GWBParser.ModiContext)):
                    recipe.append(("".join([x[0] for x in tmp]), IUPACLexer.MOD))
                elif isinstance(c, (IUPACParser.SaciContext, GWBParser.SaciContext)):
                    recipe += [(x[0], IUPACLexer.SAC) for x in tmp]
                else:
                    recipe += tmp
        if GWBLexer.AT in [x[1] for x in recipe]:
            idx = [x[1] for x in recipe].index(GWBLexer.AT)
            end = idx - 1 if idx > 0 and recipe[idx - 1][1] == GWBLexer.COMMA else idx  # check if predecessor is a comma
            recipe = recipe[:end] + recipe[idx + 2:]  # remove (,) @ NUM
        return recipe

    def add_node(self, node, parent: int = -1, config=""):
        """
        Add a new node to the network based on the name of the represented glycan.

        Args:
            node (IUPACParser.GlycanContext): Node of the parsed tree to be parsed into a monomer
            config (str): configuration to be applied to the monomer

        Returns:
            ID of the newly added node
        """
        # get the id for the node and increase the id for the next node
        node_id = self.node_id
        self.node_id += 1

        # add the node to the network and store the enum of the glycan as attribute
        recipe = self.build_recipe(node)
        # print(recipe)
        monomer, full = self.factory.create(recipe, config, tree_only=self.tree_only)
        self.g.add_node(
            node_id,
            type=monomer,
        )
        self.full &= full

        return node_id
