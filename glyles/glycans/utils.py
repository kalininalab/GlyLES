from enum import Enum
from typing import Callable

import networkx as nx
from networkx.algorithms import isomorphism
from rdkit.Chem import MolFromSmiles, FindMolChiralCenters

from itertools import permutations


class Verbosity(Enum):
    """
    define Levels of verbosity for printing different important messages
    """
    QUIET = 0
    NORMAL = 1


class Config(Enum):
    """
    Configuration enum to represent if the monomer is in alpha, beta, or undefined configuration.
    """
    UNDEF = 0
    ALPHA = 1
    BETA = 2


class Enantiomer(Enum):
    """
    Configuration of the whole monomer regarding L and D forms in terms of enantiomerism.
    """
    D = 0
    L = 1
    U = 2


class Lactole(Enum):
    """
    Specification if a monomer is a pyranose (6-ring) or a furanose (5-ring)
    """
    UNKNOWN = 0
    OPEN = 1
    FURANOSE = 5
    PYRANOSE = 6


class UnreachableError(NotImplementedError):
    """
    Represent exceptions that should arise in case a piece of code is reached that under normal circumstances should
    not be reached
    """
    pass


class ParseError(ValueError):
    """
    Represent parsing errors when reading in a glycan
    """
    pass


ketoses2 = {
    ("Ko", Lactole.PYRANOSE), ("Kde", Lactole.PYRANOSE), ("Neu", Lactole.PYRANOSE), ("Pse", Lactole.PYRANOSE),
    ("Leg", Lactole.PYRANOSE), ("Aci", Lactole.PYRANOSE), ("Kdo", Lactole.PYRANOSE), ("Dha", Lactole.PYRANOSE),
    ("Fru", Lactole.PYRANOSE), ("Sor", Lactole.PYRANOSE), ("Tag", Lactole.PYRANOSE), ("Psi", Lactole.PYRANOSE),
}


def find_rings(mol):
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) > 0:
        ring_info = [None]
        for ring in rings:
            found_ox = False
            for atom in ring:
                if mol.GetAtomWithIdx(atom).GetAtomicNum() == 8:
                    ring_info[0] = ring
                    found_ox = True
                    break
            if not found_ox:
                ring_info.append(ring)
    else:
        ring_info = rings
    return ring_info


def find_longest_c_chain(mol):
    pass


def get_rings(mol, c1_find=None):
    rings = find_rings(mol)
    if len(rings) == 0:
        if c1_find is not None:
            return c1_find(mol)
        raise ValueError("Molecule should either have a ring or define a method to find c1!")

    ox_id = -1
    for a in rings[0]:
        if mol.GetAtomWithIdx(a).GetAtomicNum() == 8:
            ox_id = a
            break
    assert ox_id != -1, "Something went wrong"

    ring = rings[0][rings[0].index(ox_id):] + rings[0][:rings[0].index(ox_id)]
    return ring


def count_chiral_matches(ring1, chiral1, ring2, chiral2):
    return sum(chiral1.get(a1, "T") == chiral2.get(a2, "T") for a1, a2 in zip(ring1, ring2))


def match_groups(mol1, a1, ring1, chiral1, mol2, a2, ring2, chiral2):
    if chiral1.get(a1, "T") != chiral2.get(a2, "T"):
        return {}

    neighbors1 = [a.GetIdx() for a in mol1.GetAtomWithIdx(a1).GetNeighbors() if a.GetIdx() not in ring1]
    neighbors2 = [a.GetIdx() for a in mol2.GetAtomWithIdx(a2).GetNeighbors() if a.GetIdx() not in ring2]

    if len(neighbors1) < len(neighbors2):
        neighbors1 += [None] * (len(neighbors2) - len(neighbors1))
    elif len(neighbors1) > len(neighbors2):
        neighbors2 += [None] * (len(neighbors1) - len(neighbors2))

    best_map = {}
    for neighbors2_perm in permutations(neighbors2):
        for n1, n2 in zip(neighbors1, neighbors2_perm):
            if mol1.GetAtomWithIdx(n1).GetAtomicNum() == mol2.GetAtomWithIdx(n2).GetAtomicNum():
                tmp_map = {n1: n2}
                tmp_map.update(match_groups(mol1, n1, ring1 + [a1], chiral1, mol2, n2, ring2 + [a2], chiral2))
                if len(tmp_map) > len(best_map):
                    best_map = tmp_map

    return best_map


def find_isomorphism(mol1: str, mol2: str):
    mol1 = MolFromSmiles(mol1)
    chiral1 = dict(FindMolChiralCenters(mol1))
    mol2 = MolFromSmiles(mol2)
    chiral2 = dict(FindMolChiralCenters(mol2))

    ring1 = get_rings(mol1)
    ring2 = get_rings(mol2)

    mapping = {}

    fwd = count_chiral_matches(ring1, chiral1, ring2, chiral2)
    bwd = count_chiral_matches(ring1, chiral1, [ring2[0]] + list(reversed(ring2[1:])), chiral2)

    print(fwd, "|", bwd)

    ring1 = list(ring1)

    if fwd >= bwd:
        for k, v in zip(ring1, ring2):
            mapping[k] = v
        ring2 = list(ring2)
    else:
        for k, v in zip(ring1, [ring2[0]] + list(reversed(ring2[1:]))):
            mapping[k] = v
        ring2 = [ring2[0]] + list(reversed(ring2[1:]))

    print(mapping)
    for a in ring1:
        mapping.update(match_groups(mol1, a, ring1, chiral1, mol2, mapping[a], ring2, chiral2))

    return mapping


def mol_to_nx(mol):
    g = nx.Graph()
    for atom in mol.GetAtoms():
        g.add_node(atom.GetIdx(), type=atom.GetAtomicNum())  # , chiral=chiral[atom.GetIdx()])
    for bond in mol.GetBonds():
        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), type=bond.GetBondTypeAsDouble())
    return g


def networkx_match_rings(mol1_ring, mol2_ring):
    matcher = isomorphism.GraphMatcher(mol1_ring, mol2_ring, node_match=lambda n1, n2: n1["type"] == n2["type"],
                                       edge_match=lambda e1, e2: e1["type"] == e2["type"])
    for mapping in matcher.subgraph_isomorphisms_iter():
        yield dict(mapping)


def networkx_fragment_isomorphism(mol1_nx, ring1, mol2_nx, ring2):
    longest_iso = {}
    mol1_no_ring = mol1_nx.subgraph(set(list(mol1_nx.nodes)).difference(set(ring1)))
    mol2_no_ring = mol2_nx.subgraph(set(list(mol2_nx.nodes)).difference(set(ring2)))
    for iso in networkx_match_rings(mol1_nx.subgraph(ring1), mol2_nx.subgraph(ring2)):
        ring_iso = {}
        for k, v in iso.items():
            k_neighbors = [x for x in mol1_nx.neighbors(k) if x not in ring1]
            v_neighbors = [x for x in mol2_nx.neighbors(v) if x not in ring2]
            if len(k_neighbors) < len(v_neighbors):
                k_neighbors += [None for _ in range((len(v_neighbors) - len(k_neighbors)))]
            elif len(v_neighbors) < len(k_neighbors):
                v_neighbors += [None for _ in range((len(k_neighbors) - len(v_neighbors)))]

            pos_iso = {}
            for v_neigh in permutations(v_neighbors):
                fg_iso = {}
                for k_n, v_n in zip(k_neighbors, v_neigh):
                    if k_n is None or v_n is None:
                        continue

                    pair_iso = {}
                    matcher = isomorphism.GraphMatcher(
                        mol1_nx.subgraph(nx.node_connected_component(mol1_no_ring, k_n).union([k])),
                        mol2_nx.subgraph(nx.node_connected_component(mol2_no_ring, v_n).union([v])),
                        node_match=lambda n1, n2: n1["type"] == n2["type"],
                        edge_match=lambda e1, e2: e1["type"] == e2["type"]
                    )
                    for mapping in matcher.subgraph_isomorphisms_iter():
                        if len(mapping) > len(pair_iso):
                            pair_iso = mapping

                    fg_iso.update(pair_iso)

                if len(fg_iso) > len(pos_iso):
                    pos_iso = fg_iso

            ring_iso.update(pos_iso)

        ring_iso.update(iso)

        if len(ring_iso) > len(longest_iso):
            longest_iso = ring_iso

    return longest_iso


def find_isomorphism_nx(mol1: str, mol2: str, c1_find: Callable = None):
    mol1_rd = MolFromSmiles(mol1)
    mol2_rd = MolFromSmiles(mol2)

    ring1 = get_rings(mol1_rd, c1_find)
    ring2 = get_rings(mol2_rd, c1_find)

    mol1_nx = mol_to_nx(mol1_rd)
    mol2_nx = mol_to_nx(mol2_rd)

    return networkx_fragment_isomorphism(mol1_nx, ring1, mol2_nx, ring2)


class Node:
    """
    Represent a node in a tree
    """

    def __init__(self, node_id, parent_id, depth, tree):
        """
        Initialize a node in of a tree
        Args:
            node_id (int): id of the node itself in the tree
            parent_id (int): id of the parent in the tree
            depth (int): depth of the node in the tree this node belongs to
            tree (Tree): reference to the tree this node belongs to
        """
        self.children = []
        self.node_id = node_id
        self.parent_id = parent_id
        self.depth = depth
        self.tree = tree

    def add_child(self, child_id):
        """
        Add a child node to this node

        Args:
            child_id (int): id of the added child node

        Returns:
            Nothing
        """
        self.children.append(child_id)

    def is_leaf(self):
        """
        Check if this node is a leaf

        Returns:
            true if this node has no children
        """
        return len(self.children) == 0

    def __str__(self):
        """
        Convert the node into a string representation

        Returns:
            String representation of this node and all its children
        """
        return "(" + str(self.node_id) + " [" + \
               ",".join([str(self.tree.nodes[child]) for child in self.children]) + "])"


class Tree:
    """
    Represent a tree with nodes
    """

    def __init__(self):
        """
        Initialize the tree as empty tree without any node
        """
        self.nodes = {}
        self.root = None

    def __str__(self):
        """
        Convert this tree into a string representation of all its nodes

        Returns:
            String representation of the tree seen from its root node
        """
        return str(self.nodes[self.root])

    def add_node(self, node_id, parent_id=-1):
        """
        Add a node to the tree

        Args:
            node_id (int): id of the new node
            parent_id (int): id of the parent of the new node

        Returns:
            Nothing
        """
        # check if the node is a root and the tree already has a root
        if parent_id == -1 and self.root is not None:
            raise ValueError("Tree cannot have two roots")
        if parent_id != -1 and parent_id not in self.nodes:
            raise ValueError("Parent node unknown")

        # add the node and eventually set the root to be that node
        if parent_id == -1:
            self.root = node_id
            self.nodes[node_id] = Node(node_id, parent_id, 0, self)
        else:
            self.nodes[node_id] = Node(node_id, parent_id, self.nodes[parent_id].depth + 1, self)

        # add the node to the list of children in its parent
        if parent_id != -1:
            self.nodes[parent_id].add_child(node_id)

    def deepest_node(self):
        """
        Find the deepest node in the tree (most distant to root)

        Returns:
            ID of the deepest node and its depth in the tree
        """
        deepest_id, deepest_depth = 0, 0
        for n_id, node in self.nodes.items():
            if node.depth > deepest_depth:
                deepest_depth = node.depth
                deepest_id = n_id
        return deepest_id, deepest_depth

    def rehang_tree(self, node_id):
        """
        Reorder the tree to start with the node with the given id in the root

        Args:
            node_id (int): is of the node to be used as new root node

        Returns:
            New tree with the specified node as root node
        """
        tree = Tree()
        stack = [(-1, node_id)]
        while len(stack) != 0:
            p_id, c_id = stack[-1]
            stack = stack[:-1]
            tree.add_node(c_id, p_id)

            children = self.nodes[c_id].children + [self.nodes[c_id].parent_id]
            for c in children:
                if c not in tree.nodes and c != -1:
                    stack.append((c_id, c))

        return tree

    def longest_chain(self, node_id=None):
        """
        Find the longest chain of nodes in the parent. This chain starts in the provided node and goes down in the
        tree. This might lead to false results of the node with the given id has two or more children.

        Args:
            node_id (int): id of the node to start in, if None, the root will be the start

        Returns:
            List of nodes along deepest path down the tree
        """
        if node_id is None:
            node_id = self.root

        if self.nodes[node_id].is_leaf():
            return [node_id]

        longest_chain = []
        for child in self.nodes[node_id].children:
            tmp = self.longest_chain(child)
            if len(tmp) > len(longest_chain):
                longest_chain = tmp

        return [node_id] + longest_chain
