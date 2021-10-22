import networkx as nx

from glyles.glycans.glycans import Glycan, Chirality


class SMILES:
    def __init__(self):
        pass

    def write(self, g, source):
        # TODO: Reimplement DFS to get a search along the inner ring in clockwise manner
        dfs_tree = nx.algorithms.dfs_tree(g, source)
        path = nx.algorithms.dag.dag_longest_path(dfs_tree)
        return self.__construct(g, dfs_tree, path, source, first=True)

    def read(self):
        pass

    def __construct(self, g, tree, path, node, first=False):
        """

        Args:
            g (nx.Graph): graph of the molecule
            tree (nx.DiGraph): dfs-tree on that molecule
            path (list): longest path in tree
            node (int): node id the construction is currently in

        Returns:
            SMILES representation of the subtree of the current node in tree
        """
        if g.nodes[node]["chiral"] != Chirality.N:  # chirality
            output = "[C{value}H]".format(value=g.nodes[node]["chiral"].value)
        else:
            output = g.nodes[node]["type"].value

        if (len(tree.edges(node)) + 1) != len(g.edges(node)) or first:  # and node in path:
            output += "1"

        children = list(tree.edges(node))
        if len(children) == 0:  # leaf
            return output
        else:  # children
            path_one = -1
            for c in [x[1] for x in children]:
                if c in path:
                    path_one = c
                    continue
                output += "({smiles})".format(smiles=self.__construct(g, tree, path, c))
            return output + ("" if path_one == -1 else self.__construct(g, tree, path, path_one))


print(SMILES().write(Glycan.GLC.structure(), 1))
