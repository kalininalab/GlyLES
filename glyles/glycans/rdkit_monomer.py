import numpy as np
from rdkit.Chem import MolFromSmiles, MolToSmiles, GetAdjacencyMatrix, GetShortestPath

from glyles.glycans.monomer import Monomer


class RDKitMonomer(Monomer):
    def __init__(self, origin=None, **kwargs):
        """
        Initialize the monomer using the super method. Additionally, some fields are initialized to describe the
        structure of the monomer according to the specification of the monomer-parent class

        Args:
            origin (Monomer): Other monomer to use to initialize this object
            **kwargs: arguments to initialize monomer if object is None. Must include name, SMILES, and config
        """
        super(RDKitMonomer, self).__init__(origin, **kwargs)
        if isinstance(origin, RDKitMonomer):
            self._adjacency = origin.get_adjacency()
            self._ring_info = origin.get_ring_info()
            self._x = origin.get_features()
        else:
            self._adjacency = None
            self._ring_info = None
            self._x = None
            self.__get_structure()

    def get_adjacency(self):
        """
        Get the adjacency matrix of the atoms in this monomer.

        Returns:
            Adjacency matrix of all non-hydrogen atoms in this monomer
        """
        return self._adjacency

    def get_ring_info(self):
        """
        Get information of all atom-ids in the rings of this monomer.

        Returns:
            Tuple of tuples with the atom-ids from rdkit in the monomer
        """
        return self._ring_info

    def get_features(self):
        """
        Get a feature matrix from this monomer. The features contain information about atom type, ids according to the
        specification in monomer.py, and ring memberships.

        Returns:
            A numpy array of shape Nx3 containing the extracted features for all atoms in this molecule
        """
        return self._x

    def get_dummy_atoms(self):
        """
        Specify some dummy atoms that are used to mark oxygen atoms that will participate in bindings between glycans.
        Here, the atoms will be replaced by instances of the atom enum that are used to define the type of the atom in
        the nodes of the networkx representation of the monomer molecules.

        Returns:
            Two lists:
                * one with atoms that are given as the "atom" argument to the mark-method to replace the oxygen atoms
                * the string representation of the atoms from above, i.e. how the atoms above will be represented in a
                  SMILES string
        """
        return [34, 52, 84], ["[SeH]", "[TeH]", "[PoH]"]

    def root_atom_id(self, binding_c_id):
        """
        Get ID of atom that will bind the parent monomer in the glycan. This ID will be given as root argument to
        the to_smiles method.

        Args:
            binding_c_id (int): Integer at which c-position this monomer binds its parent

        Returns:
            id of the atom that binds to the parent, -1 if the root cannot be found
        """
        return self.__find_oxygen(binding_c_id)

    def mark(self, position, atom):
        """
        Mark the oxygen atom linked to the carbon atom at the given position ready to participate in the bounding.
        Marking here works based on replacing the oxygen-group bound to the carbon atom at the given position with the
        atom enum instance also provided in the arguments.

        Args:
            position (int): id of the carbon atom whose oxygen atom will from the binding
            atom (object): atom to replace the binding oxygen with

        Returns:
            Nothing
        """
        idx = self.__find_oxygen(position)
        self.__get_structure().GetAtomWithIdx(idx).SetAtomicNum(atom)
        self._x[idx, 0] = atom

    def to_smiles(self, root, ring_index):
        """
        Convert this monomer into a SMILES string representation.
        Use the implementation of the SMILES algorithm fitted to the needs of glycans.

        Args:
            root (int): index of the root atom
            ring_index (int): index of the rings in the atom

        Returns:
            SMILES string representation of this molecule
        """
        smiles = MolToSmiles(self.__get_structure(), rootedAtAtom=root)
        return "".join([(str(int(c) + ring_index) if c.isdigit() else c) for c in smiles])

    def __get_structure(self):
        """
        Compute and save the structure of this glycan.

        Returns:
            rdkit molecule representing the structure of the glycan as a graph of its non-hydrogen atoms.
        """
        if self._structure is None:
            # read the structure from the SMILES string
            self._structure = MolFromSmiles(self._smiles)

            # extract some further information from the molecule to not operate always on the molecule
            self._adjacency = GetAdjacencyMatrix(self._structure)
            self._ring_info = self._structure.GetRingInfo().AtomRings()
            self._x = np.zeros((self._adjacency.shape[0], 3))

            c_atoms, ringo = [], -1
            # extract some information form the molecule
            for i in range(self._adjacency.shape[0]):
                atom = self._structure.GetAtomWithIdx(i)

                # store the atom type
                self._x[i, 0] = atom.GetAtomicNum()
                if self._x[i, 0] == 6 and i in self._ring_info[0]:
                    c_atoms.append(int(i))

                # if the atom is part of any ring, store the number of that ring
                for r in range(len(self._ring_info)):
                    if i in self._ring_info[r]:
                        self._x[i, 2] = r + 1

                # identify the oxygen atom in the main ring and set its id to 10
                if self._x[i, 2] == 1 and self._x[i, 0] == 8:
                    self._x[i, 1] = 10
                    ringo = i

            self.__enumerate_c_atoms2(c_atoms, ringo)
        return self._structure

    class Tree:
        class Node:
            def __init__(self, node_id, parent_id, depth, tree):
                self.children = []
                self.node_id = node_id
                self.parent_id = parent_id
                self.depth = depth
                self.tree = tree

            def add_child(self, child_id):
                self.children.append(child_id)

            def is_leaf(self):
                return len(self.children) == 0

            def __str__(self):
                return "(" + str(self.node_id) + " [" + ",".join([str(self.tree.nodes[child]) for child in self.children]) + "])"

        def __init__(self):
            self.nodes = {}
            self.root = None

        def __str__(self):
            return str(self.nodes[self.root])

        def add_node(self, node_id, parent_id=-1):
            if parent_id == -1 and self.root is not None:
                raise ValueError("Tree cannot have two roots")
            if parent_id == -1:
                self.root = node_id
                self.nodes[node_id] = RDKitMonomer.Tree.Node(node_id, parent_id, 0, self)
            else:
                self.nodes[node_id] = RDKitMonomer.Tree.Node(node_id, parent_id, self.nodes[parent_id].depth + 1, self)
            if parent_id != -1:
                self.nodes[parent_id].add_child(node_id)

        def deepest_node(self):
            deepest_id, deepest_depth = 0, 0
            for n_id, node in self.nodes.items():
                if node.depth > deepest_depth:
                    deepest_depth = node.depth
                    deepest_id = n_id
            return deepest_id, deepest_depth

        def rehang_tree(self, node_id):
            tree = RDKitMonomer.Tree()
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

    def __equidistant(self, start, end, ringo):
        pass

    def __evaluate_distance(self, start, end, ringo):
        adj = self._adjacency.copy()

        while adj[start, ringo] == 0 and adj[end, ringo] == 0:
            adj *= self._adjacency

        if adj[start, ringo] > 0 and adj[end, ringo] > 0:
            return self.__equidistant(start, end, ringo)
        elif adj[start, ringo] > 0:
            return True
        return False

    def __enumerate_c_atoms2(self, c_atoms, ringo):
        # find longest c-chain
        c_tree = RDKitMonomer.Tree()
        stack = [(-1, c_atoms[0])]
        while len(stack) != 0:
            p_id, c_id = stack[-1]
            stack = stack[:-1]
            c_tree.add_node(c_id, p_id)

            children = np.argwhere(self._adjacency[c_id] & (self._x[:, 0] == 6))
            for c in children:
                if int(c) not in c_tree.nodes:
                    stack.append((c_id, int(c)))

        deepest_id, _ = c_tree.deepest_node()
        c_tree = c_tree.rehang_tree(deepest_id)
        longest_c_chain = c_tree.longest_chain()

        # find ends
        start, end = longest_c_chain[0], longest_c_chain[-1]

        # check conditions
        start_o_conn = any(self._x[np.where(self._adjacency[start, :]), 0].flatten() == 6)
        end_o_conn = any(self._x[np.where(self._adjacency[end, :]), 0].flatten() == 6)

        # decide on c1
        if start_o_conn and end_o_conn:
            if not self.__evaluate_distance(start, end, ringo):
                longest_c_chain = reversed(longest_c_chain)
        elif end_o_conn:
            longest_c_chain = reversed(longest_c_chain)

        # enumerate along chain
        c_count = 0
        for c in longest_c_chain:
            c_count += 1
            self._x[c, 1] = c_count

    def __find_oxygen(self, binding_c_id):
        """
        Find the oxygen atom that binds to the carbon atom with the provided id. The returned id may not refer to an
        oxygen atom in the ring of the monomer as this cannot bind anything.

        Args:
            binding_c_id (int): id of the carbon atom that participates in a binding, and we need to find the oxygen
                from

        Returns:
            id referring to the oxygen binding the provided carbon atom and may participate in a glycan-binding.
        """
        # first find the rdkit id of the carbon atom that should bind to something
        position = np.argwhere(self._x[:, 1] == binding_c_id).squeeze()

        # then find the candidates. There should be exactly one element in the resulting array
        candidates = np.argwhere((self._adjacency[position, :] == 1) &
                                 (self._x[:, 0] == 8) & (self._x[:, 2] != 1)).squeeze()
        return int(candidates)

    @staticmethod
    def from_string(mono):
        """
        Convert the string representation of the monomer into an instance.

        Args:
            mono (str): string encoding of a monomer

        Returns:
            RDKit instance according to the string
        """
        tmp = Monomer.from_string(mono)
        return RDKitMonomer(origin=tmp)
