from enum import Enum


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
