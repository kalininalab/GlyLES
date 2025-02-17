import logging

from antlr4 import CommonTokenStream, InputStream
import networkx as nx

from glyles.glycans.poly.glycan import Glycan
from glyles.glycans.poly.gwb_walker import GWBTreeWalker
from glyles.glycans.poly.merger import Merger
from glyles.glycans.utils import ParseError
from glyles.gwb.GWBLexer import GWBLexer
from glyles.gwb.GWBParser import GWBParser


def graph_to_string_int(graph: nx.DiGraph, node2label: callable) -> str:
    "Convert glycan graph back to IUPAC-condensed format"
    def assign_depth(G: nx.DiGraph, idx: int) -> float:
        "Assign depth to each node in the graph recursively"
        min_depth = float("inf")
        children = list(G.neighbors(idx))
        if len(children) == 0:  # if it's a leaf node, the depth is zero
            G.nodes[idx]["depth"] = 0
            return 0
        for child in children:  # if it's not a leaf node, the depth is the minimum depth of its children + 1
            min_depth = min(assign_depth(G, child) + 1, min_depth)
        G.nodes[idx]["depth"] = min_depth
        return min_depth

    def dfs_to_string(G: nx.DiGraph, idx: int, parent: int) -> str:
        "Convert the DFS tree of a glycan graph to IUPAC-condensed format recursively"
        output = node2label(graph.nodes[idx]["type"])
        if parent != -1:  # put the linkage information in the output string
            output += str(graph.get_edge_data(parent, idx)["type"])
        elif "float_edge" in graph.nodes[idx]:
            output += str(graph.nodes[idx]["float_edge"])
        # sort kids from shallow to deep to have deepest as main branch in the end
        children = sorted([(dfs_to_string(G, child, idx), G.nodes[child]["depth"]) for child in G.neighbors(idx)], key=lambda x: (x[1], x[0]))
        if len(children) == 0:
            return output
        for child_iupac, _ in children[:-1]:  # iterate over all children except the last one and put them in branching brackets
            output = "[" + child_iupac + "]" + output
        # put the last child in front of the output, without brackets
        output = children[-1][0] + output
        return output

    # get the root node index, assuming the root node has the lowest (highest) index
    root_idx = min(list(graph.nodes.keys()))  # max(list(graph.nodes.keys()))
    # get the DFS tree of the graph
    dfs = nx.dfs_tree(graph, root_idx)
    # assign depth to each node
    assign_depth(dfs, root_idx)
    # convert the DFS tree to IUPAC-condensed format
    return dfs_to_string(dfs, root_idx, -1)


def graph_to_string(graph: nx.Graph, node2label: callable) -> str: # IUPAC-condensed glycan string
    "Convert glycan graph back to IUPAC-condensed format, handling disconnected components"
    if nx.number_connected_components(graph.to_undirected()) > 1:
        parts = [graph.subgraph(sorted(c)) for c in nx.connected_components(graph.to_undirected())]
        main_iupac = ""
        floats = []
        for p in parts:
            iupac = graph_to_string_int(p, node2label)
            if iupac[-1] == ")":
                floats.append(iupac)
            else:
                main_iupac = iupac
        return "{" + "}{".join(sorted(floats, key=lambda x: (len(x), x))) + "}" + main_iupac
    else:
        return graph_to_string_int(graph, node2label)


class GWBGlycan(Glycan):
    def to_iupac(self, slim: bool = False) -> str:
        return graph_to_string(self.parse_tree, lambda x: x.get_name(mode="slim"))
    
    def parse(self):
        try:
            if "$" in self.iupac:
                self.iupac = self.iupac[:self.iupac.index("$")]
            self.iupac += "$"
            stream = InputStream(data=self.iupac)
            lexer = GWBLexer(stream)
            token = CommonTokenStream(lexer)
            parser = GWBParser(token)

            try:
                self.grammar_tree = parser.start()
            except ParseError as pe:
                self.parse_tree = None
                self.glycan_smiles = ""
                raise pe from None
            self.parse_tree, self.tree_full = GWBTreeWalker(self.factory, self.tree_only).parse(self.grammar_tree)

            # if the glycan should be parsed immediately, do so
            if not self.tree_only and self.tree_full and self.full:
                self.glycan_smiles = Merger(self.factory).merge(self.parse_tree, self.root_orientation, start=self.start)
                # catch any exception at glycan level to not destroy the whole pipeline because of one mis-formed glycan
        except ParseError as e:
            msg = e.__str__().replace("\n", " ")
            logging.error(f"A parsing error occurred with \"{self.iupac}\": Error message: {msg}")
            self.glycan_smiles = ""
        except Exception as e:
            msg = e.__str__().replace("\n", " ")
            logging.error(f"An unexpected exception occurred with \"{self.iupac}\". This glycan cannot be parsed. "
                          f"Error message: {msg}")
            self.glycan_smiles = ""
            raise e


if __name__ == "__main__":
    print(GWBGlycan("redEnd--??1D-Glc,p--4b1D-Gal,p(--3a2D-NeuGc,p@270)--4b1D-GalNAc,p$").to_iupac(slim=True))
