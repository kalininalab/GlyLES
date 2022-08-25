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
    """
    Determine the rings of a molecule in RDKit.

    Args:
        mol (rdkit.Molecule): molecule to determine the rings for

    Returns:
        List (possibly empty) of rings of the input molecule. The list will contain a list of RDKit IDs of every atom
        in the respective rings
    """
    # get the rings from the molecule
    rings = mol.GetRingInfo().AtomRings()

    if len(rings) > 0:

        # if there are rings, reorder them to have the main-ring (with the oxygen) at first position
        ring_info = [None]
        for ring in rings:
            found_ox = False

            # if the ring contains an oxygen atom, it's the first/main ring
            for atom in ring:
                if mol.GetAtomWithIdx(atom).GetAtomicNum() == 8:
                    ring_info[0] = ring
                    found_ox = True
                    break

            # otherwise, just append it to the list of rings
            if not found_ox:
                ring_info.append(ring)

    # if there are no rings, return an empty list
    else:
        ring_info = rings

    return ring_info


def get_rings(mol, name, c1_find=None):
    """
    Get the longest carbon chain in the molecule starting from the C1 atom based on the rings of the molecule.

    Args:
        mol (rdkit.Molecule): molecule to determine the carbon list for
        name (str): Name of the monosaccharide we're working with
        c1_find (Callable): Optional function to determine C1 for non-ring monosaccharides

    Returns:
        List of the longest carbon chain starting with C1
    """
    # find the rings
    rings = find_rings(mol)
    if len(rings) == 0 or rings[0] is None or name == "Inositol":
        if c1_find is not None:
            return c1_find(mol)
        raise ValueError("Molecule should either have a ring or define a method to find c1!")

    # identify the oxygen in the first/main ring
    ox_id = -1
    for a in rings[0]:
        if mol.GetAtomWithIdx(a).GetAtomicNum() == 8:
            ox_id = a
            break
    if ox_id == -1:
        raise UnreachableError("First/Main should always have an oxygen atom in it!")

    # reorient the carbon list to start with the oxygen
    ring = rings[0][rings[0].index(ox_id):] + rings[0][:rings[0].index(ox_id)]

    return ring


def mol_to_nx(mol):
    """
    Convert a molecule to a networkx graph with atom types and bond types as only features.

    Args:
        mol (rdkit.Molecule): molecule to be converted

    Returns:
        networkx graph representing a molecule
    """
    g = nx.Graph()

    # add all the nodes/atoms
    for atom in mol.GetAtoms():
        g.add_node(atom.GetIdx(), type=atom.GetAtomicNum())

    # add all the bonds/edges
    for bond in mol.GetBonds():
        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), type=bond.GetBondTypeAsDouble())
    return g


def networkx_match_rings(mol1_ring, mol2_ring):
    """
    Map the rings onto each other and yield every possible isomorphism between the rings

    Args:
        mol1_ring (networkx.Graph): graph representing the fist molecules ring
        mol2_ring (networkx.Graph): graph representing the second molecules ring

    Returns:
        Mappings between both rings corresponding to isomorphisms between both rings
    """
    # define a matcher for both graphs with the requirement matching all bond types and edge types for an isomorphism
    matcher = isomorphism.GraphMatcher(mol1_ring, mol2_ring, node_match=lambda n1, n2: n1["type"] == n2["type"],
                                       edge_match=lambda e1, e2: e1["type"] == e2["type"])

    # iterate over all possible isomorphism fulfilling previously mentioned requirement and yield them
    for mapping in matcher.subgraph_isomorphisms_iter():
        yield dict(mapping)


def networkx_fragment_isomorphism(mol1_nx, ring1, mol2_nx, ring2):
    """
    Find the biggest isomorphism between the two molecules represented in the graphs.

    Args:
        mol1_nx (networkx.Graph): graph of molecule1
        ring1 (List[int]): main-ring of molecule 1
        mol2_nx (networkx.Graph): graph of molecule 2
        ring2 (List[int]): ring of molecule 2

    Returns:
        Mapping from mol1's RDKIT IDs to mol2's RDKit IDs
    """
    # keep track of the biggest isomorphism between both molecules
    longest_iso = {}

    # create the graphs of the side chains (it's actually a collection of trees)
    mol1_no_ring = mol1_nx.subgraph(set(list(mol1_nx.nodes)).difference(set(ring1)))
    mol2_no_ring = mol2_nx.subgraph(set(list(mol2_nx.nodes)).difference(set(ring2)))

    # iterate over all isomorphisms between just the two rings
    for iso in networkx_match_rings(mol1_nx.subgraph(ring1), mol2_nx.subgraph(ring2)):
        ring_iso = {}

        # iterate over every matched position in the isomorphism
        for k, v in iso.items():

            # extract all neighbors at the given positions in the molecules
            k_neighbors = [x for x in mol1_nx.neighbors(k) if x not in ring1]
            v_neighbors = [x for x in mol2_nx.neighbors(v) if x not in ring2]

            # pad the neighborhood lists with Nones to have same length
            if len(k_neighbors) < len(v_neighbors):
                k_neighbors += [None for _ in range((len(v_neighbors) - len(k_neighbors)))]
            elif len(v_neighbors) < len(k_neighbors):
                v_neighbors += [None for _ in range((len(k_neighbors) - len(v_neighbors)))]

            pos_iso = {}

            # iterate over all permutation of the latter list...
            for v_neigh in permutations(v_neighbors):
                fg_iso = {}

                # ... to find the best isomorphism between the functional groups of Ck from mol1 and Cv from mol2
                for k_n, v_n in zip(k_neighbors, v_neigh):

                    # if one of the two side-chains is None, there is no isomorphism, so continue to next match
                    if k_n is None or v_n is None:
                        continue

                    pair_iso = {}

                    # define a matcher between the side chain trees ...
                    matcher = isomorphism.GraphMatcher(
                        mol1_nx.subgraph(nx.node_connected_component(mol1_no_ring, k_n).union([k])),
                        mol2_nx.subgraph(nx.node_connected_component(mol2_no_ring, v_n).union([v])),
                        node_match=lambda n1, n2: n1["type"] == n2["type"],
                        edge_match=lambda e1, e2: e1["type"] == e2["type"]
                    )

                    # ... and extract their best match
                    for mapping in matcher.subgraph_isomorphisms_iter():
                        if len(mapping) > len(pair_iso):
                            pair_iso = mapping

                    # put the best match of the side chains into the matching for the position isomorphisms
                    fg_iso.update(pair_iso)

                # keep track of the best isomorphisms between (potentially) different side chains at the same position
                if len(fg_iso) > len(pos_iso):
                    pos_iso = fg_iso

            # for each position, save the best isomorphism between the side chain trees
            ring_iso.update(pos_iso)

        # add the isomorphism of the side chains with the isomorphism of the ring
        ring_iso.update(iso)

        # save the best/biggest isomorphism between the two molecules
        if len(ring_iso) > len(longest_iso):
            longest_iso = ring_iso

    return longest_iso


def find_isomorphism_nx(mol1, mol2, name, c1_find=None):
    """
    Find an isomorphism between two molecules. The problem is an NP problem, normally. Here, we have some domain
    knowledge based on how and when this method is called.
    First, there is an isomorphism of at least the two rings. Second, the isomorphism can be extended by the attached
    functional groups. These functional groups are tree-like structured and not intersecting which makes it much easier
    to find isomorphisms between them.

    Args:
        mol1 (str): SMILES of the monosaccharide with functional groups
        mol2 (str): SMILES of the root monosaccharide
        name (str): Name of the monosaccharide we're working with
        c1_find (Callable): Optional method to provide to find C1 in special molecules

    Returns:
        Mapping from mol1's RDKIT IDs to mol2's RDKit IDs
    """
    # generate the RDKit molecules from the SMILES strings
    mol1_rd = MolFromSmiles(mol1)
    mol2_rd = MolFromSmiles(mol2)

    # identify the rings or the longest carbon-cain starting with
    ring1 = get_rings(mol1_rd, name, c1_find)
    ring2 = get_rings(mol2_rd, name, c1_find)

    # convert the molecules to networkx graphs for easier accessibility
    mol1_nx = mol_to_nx(mol1_rd)
    mol2_nx = mol_to_nx(mol2_rd)

    # actually compute the isomorphism between the two molecules and return it
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
