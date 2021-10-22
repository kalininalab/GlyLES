import networkx as nx

from glyles.glycans.glycans import Glycan, Chirality


class DFS:
    def dfs_tree(self, g, source):
        dfst = nx.DiGraph()
        self.__dfs(g, dfst, None, source)

        return dfst

    def __dfs(self, g, dfst, parent, node):
        dfst.add_node(node)

        if parent is not None:
            dfst.add_edge(parent, node)

        children = sorted([x[1] for x in list(g.edges(node))], key=lambda x: (not g.nodes[x]["ring"], x))
        for child in children:
            if child not in dfst.nodes:
                self.__dfs(g, dfst, node, child)


class SMILES:
    def __init__(self):
        pass

    def write(self, g, source):
        dfs_tree = DFS().dfs_tree(g, source)
        return self.__construct(g, dfs_tree, source, first=True)

    def read(self):
        pass

    def __construct(self, g, tree, node, first=False):
        """

        Args:
            g (nx.Graph): graph of the molecule
            tree (nx.DiGraph): dfs-tree on that molecule
            node (int): node id the construction is currently in

        Returns:
            SMILES representation of the subtree of the current node in tree
        """
        if g.nodes[node]["chiral"] != Chirality.NONE:  # chirality
            output = "[C{value}H]".format(value=("@@" if g.nodes[node]["chiral"] == Chirality.DOWN else "@"))
        else:
            output = g.nodes[node]["type"].value

        if (len(tree.edges(node)) + 1) != len(g.edges(node)) or first:  # and node in path:
            output += "1"

        children = list(tree.edges(node))
        if len(children) == 0:  # leaf
            return output
        else:  # children
            for c in [x[1] for x in children]:
                output += "({smiles})".format(smiles=self.__construct(g, tree, c))
            return output


print(SMILES().write(Glycan.GLC.structure(), 1))
