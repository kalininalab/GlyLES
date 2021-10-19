import networkx as nx
from antlr4 import *
from antlr4.tree.Tree import ParseTree, TerminalNodeImpl

from glyles.glycans import glycans
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.grammar.GlycanListener import GlycanListener
from glyles.grammar.GlycanParser import GlycanParser


class TreeWalker:
    def __init__(self):
        self.g = nx.Graph()
        self.node_id = 0

    def parse(self, t: ParseTree, init):
        node_id = self.__add_node(init)
        self.__walk(next(t.getChildren()), node_id)
        return self.g

    def __walk(self, t: ParseTree, parent):
        if isinstance(t, ErrorNode):
            raise NotImplementedError("ErrorNodes not implemented yet!")
        elif isinstance(t, TerminalNode):
            raise NotImplementedError("Terminal nodes should be unreachable!")

        children = list(t.getChildren())

        if len(children) == 2:  # [sac con]
            node_id = self.__add_node(children[0].symbol.text)
            self.__add_edge(parent, node_id, children[1])
            return node_id
        elif len(children) == 3 and isinstance(children[2], GlycanParser.BranchContext):  # [sac con (sac con ...)]
            parent = self.__walk(children[2], parent)
            node_id = self.__add_node(children[0].symbol.text)
            self.__add_edge(parent, node_id, children[1])
            return node_id
        elif len(children) == 3 and isinstance(children[1], GlycanParser.BranchContext):  # [ [ branch ] ]
            self.__walk(children[1], parent)
            return parent
        elif len(children) == 7:
            node_id = self.__add_node(children[5].symbol.text)
            self.__add_edge(parent, node_id, children[6])
            self.__walk(children[3], node_id)
            node_id2 = self.__add_node(children[0].symbol.text)
            self.__add_edge(node_id, node_id2, children[1])
            return node_id2
        raise NotImplementedError("This should be unreachable")

    def __add_edge(self, parent, child, con):
        # children = list(con.getChildren())
        # self.g.add_edge(parent, child, attr={(children[1].symbol.text, children[2].symbol.text, children[4].symbol.text)})
        self.g.add_edge(parent, child, attr={(con.symbol.text)})

    def __add_node(self, name):
        node_id = self.node_id
        self.g.add_node(node_id, attr={"type": glycans.from_string(name)})
        self.node_id += 1
        return node_id


def parse(smiles: str):
    if "(" not in smiles:
        raise NotImplementedError("Monosaccharides are not implemented yet!")

    bracket_index = max(
        smiles.rindex(")") if ")" in smiles else -1,
        smiles.rindex("]") if "]" in smiles else -1,
    ) + 1
    init = smiles[bracket_index:]
    smiles = smiles[:bracket_index]

    stream = InputStream(data=smiles)

    lexer = GlycanLexer(stream)
    token = CommonTokenStream(lexer)
    parser = GlycanParser(token)
    tree = parser.start()

    if True:
        walker = TreeWalker()
        g = walker.parse(tree, init)
        nx.drawing.nx_pydot.write_dot(g, f"./{smiles}{init}.dot")
