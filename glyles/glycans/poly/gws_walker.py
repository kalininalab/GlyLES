import networkx as nx
from antlr4 import ErrorNode, TerminalNode

from glyles.glycans.poly.walker import TreeWalker
from glyles.glycans.utils import UnreachableError
from glyles.gws.GWSLexer import GWSLexer
from glyles.gws.GWSParser import GWSParser


if not hasattr(GWSLexer, "MOD"):
    GWSLexer.MOD = GWSLexer.REDEND + 1


class GWSTreeWalker(TreeWalker):
    def __init__(self, factory, tree_only):
        super().__init__(factory, tree_only)
    
    def parse(self, t):
        for c in filter(lambda x: not isinstance(x, (TerminalNode, ErrorNode)), t.getChildren()):
            if isinstance(c, GWSParser.BranchContext):
                self.walk(c, -1)
            elif isinstance(c, GWSParser.FloatContext):
                self.walk(c, -1)
        return self.g, self.full and len(nx.connected_components(self.g.to_undirected())) == 1
    
    def walk(self, t, parent):
        """
        Args:
            t (GWSParser.BranchContext): The branch to parse.
        """
        if isinstance(t, ErrorNode):
            raise UnreachableError("ErrorNodes in parsing!")
        elif isinstance(t, TerminalNode):
            raise UnreachableError("TerminalNodes in parsing!")
        
        children = list(t.getChildren())
        if len(children) == 2:  # con deriv
            node_id = self.add_node(children[1])
            self.full &= self.add_edge(parent, node_id, children[0])
            return node_id
        elif len(children) == 3 and isinstance(children[1], GWSParser.BranchContext):  # ( branch )
            self.walk(children[1], parent)
            return parent
        elif len(children) % 3 == 0:  # con deriv (LPAR branch RPAR)* branch
            node_id = self.add_node(children[1])
            if len(children) > 3:
                for i in range(3, len(children), 3):
                    self.walk(children[i], node_id)  # link branches to node_id, which is the upstream monosaccharide
            self.walk(children[-1], node_id)
            self.full &= self.add_edge(parent, node_id, children[0])
            return node_id
        raise UnreachableError("Unexpected number of children in branch!")

    def add_edge(self, parent, child, con) -> bool:
        if parent == -1:
            return True
        cs = [self.context2str(c) for c in list(con.getChildren())[1:]][::-1]
        con = f"({cs[0]}{cs[1]}-{cs[2]})"
        self.g.add_edge(parent, child, type=con)
        return "?" not in con and "/" not in con
    
    def add_node(self, node, config=""):
        """
        Add a new node to the network based on the name of the represented glycan.

        Args:
            node (GlycanParser.GlycanContext): Node of the parsed tree to be parsed into a monomer
            config (str): configuration to be applied to the monomer

        Returns:
            ID of the newly added node
        """
        # get the id for the node and increase the id for the next node
        node_id = self.node_id
        self.node_id += 1

        # add the node to the network and store the enum of the glycan as attribute
        recipe = self.build_recipe(node)
        recipe = [recipe[0]] + recipe[1:-2] + [recipe[-1]]
        # print(recipe)
        monomer, full = self.factory.create(recipe, config, tree_only=self.tree_only)
        self.g.add_node(
            node_id,
            type=monomer,
        )
        self.full &= full

        return node_id
