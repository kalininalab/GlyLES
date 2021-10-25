import networkx as nx

from glyles.glycans.glycans import Glycan, Chirality, Atom
from glyles.grammar.parse import parse


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
        TODO: Test from non-circular starting points
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


class Merger:
    def merge(self, t):
        self.__mark(t, 0)
        return self.__merge(t, 0, 1)

    def __mark(self, t, node):
        children = [x[1] for x in t.edges(node)]

        if len(children) == 0:  # leaf
            return
        if len(children) > 2:
            raise NotImplementedError("Glycans with maximal branching factor not implemented.")

        for child, atom in zip(children, [Atom.Y, Atom.Z]):
            binding = list(t.get_edge_data(node, child)["type"])
            # t.nodes[child].structure().nodes[binding[2]]  # davon die OH gruppe als start für SMILES!
            # t.nodes[node].structure().nodes[binding[4]]  # davon die OH gruppe duch X ersetzen
            monomer = t.nodes[node]["type"].structure()
            for _, x in monomer.edges(int(binding[4])):   # davon die OH gruppe duch X ersetzen
                if monomer.nodes[x]["type"] == Atom.O and 11 <= x <= 15:
                    monomer.nodes[x]["type"] = atom
                    break

    def __merge(self, t, node, start):
        children = [x[1] for x in t.edges(node)]
        me = SMILES().write(t.nodes[node]["type"].structure(), start)

        if len(children) == 0:  # leaf
            return me
        if len(children) > 2:
            raise NotImplementedError("Glycans with maximal branching factor not implemented.")

        for child, atom in zip(children, ["Y", "Z"]):
            binding = list(t.get_edge_data(node, child)["type"])
            monomer = t.nodes[node]["type"].structure()
            child_start = -1

            for _, x in monomer.edges(int(binding[2])):  # davon die OH gruppe duch X ersetzen
                if monomer.nodes[x]["type"] == Atom.O and 11 <= x <= 15:
                    child_start = x
                    break

            if child_start == -1:
                raise ValueError("SMILES cannot computed from atom -1!")

            child_smiles = self.__merge(t, child, child_start)
            me = me.replace(atom, child_smiles)

        return me


# print(SMILES().write(Glycan.GLC.structure(), 1))
g = parse("Man(a1-4)[Fru(a1-3)]Gal")
s = Merger().merge(g)
print(s)
