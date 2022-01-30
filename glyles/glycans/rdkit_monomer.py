import warnings

import numpy as np
from rdkit.Chem import MolFromSmiles, MolToSmiles, GetAdjacencyMatrix
from rdkit.Chem.rdchem import Atom, EditableMol, BondType

from glyles.glycans.monomer import Monomer
from glyles.glycans.utils import UnreachableError, Tree, ModificationNotImplementedWarning
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.glycans.rdkit_reactor import Reactor

'''def not_implemented_message(mod):
    print(f"ModificationNotImplementedWarning: {mod} Modification not implemented. The returned molecule will not have "
          f"this modification")'''


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
            self.adjacency = origin.get_adjacency()
            self.ring_info = origin.get_ring_info()
            self.x = origin.get_features()
        else:
            self.adjacency = None
            self.ring_info = None
            self.x = None
            self.get_structure()

    def alpha(self, factory):
        """
        Return this monosaccharide in its alpha conformation.

        Args:
            factory (MonomerFactory): factory instance to be used to generate the monomers

        Returns:
            Monomer in alpha conformation
        """
        recipe = [(v, t) for v, t in self.recipe if t != GlycanLexer.TYPE]
        recipe.append(('a', GlycanLexer.TYPE))
        return RDKitMonomer(factory.create(recipe))

    def beta(self, factory):
        """
        Return this monosaccharide in its beta conformation.

        Args:
            factory (MonomerFactory): factory instance to be used to generate the monomers

        Returns:
            Monomer in beta conformation
        """
        recipe = [(v, t) for v, t in self.recipe if t != GlycanLexer.TYPE]
        recipe.append(('b', GlycanLexer.TYPE))
        return RDKitMonomer(factory.create(recipe))

    def undefined(self, factory):
        """
        Return this monosaccharide in undefined conformation, the first carbon ring-atom will have unspecified.
        chirality.

        Args:
            factory (MonomerFactory): factory instance to be used to generate the monomers

        Returns:
            Monomer in undefined conformation
        """
        recipe = [(v, t) for v, t in self.recipe if t != GlycanLexer.TYPE]
        return RDKitMonomer(factory.create(recipe))

    def get_adjacency(self):
        """
        Get the adjacency matrix of the atoms in this monomer.

        Returns:
            Adjacency matrix of all non-hydrogen atoms in this monomer
        """
        if self.adjacency is None:
            self.get_structure()
        return self.adjacency

    def get_ring_info(self):
        """
        Get information of all atom-ids in the rings of this monomer.

        Returns:
            Tuple of tuples with the atom-ids from rdkit in the monomer
        """
        if self.ring_info is None:
            self.get_structure()
        return self.ring_info

    def get_features(self):
        """
        Get a feature matrix from this monomer. The features contain information about atom type, ids according to the
        specification in monomer.py, and ring memberships.

        Returns:
            A numpy array of shape Nx3 containing the extracted features for all atoms in this molecule
        """
        if self.x is None:
            self.get_structure()
        return self.x

    def get_dummy_atoms(self):
        """
        Specify some dummy atoms that are used to mark oxygen atoms that will participate in bindings between glycans.
        Here, the atoms will be replaced by instances of the atom enum that are used to define the type of the atom in
        the nodes of the networkx representation of the monomer molecules.
        TODO: Extend this to N- and C-glycosidic bonds

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
        return self.find_oxygen(binding_c_id)

    def mark(self, position, atom):
        """
        Mark the oxygen atom linked to the carbon atom at the given position ready to participate in the bounding.
        Marking here works based on replacing the oxygen-group bound to the carbon atom at the given position with the
        atom enum instance also provided in the arguments.

        Args:
            position (int): id of the carbon atom whose oxygen atom will from the binding
            atom (int): atom to replace the binding oxygen with

        Returns:
            Nothing
        """
        idx = self.find_oxygen(position)
        self.get_structure().GetAtomWithIdx(idx).SetAtomicNum(atom)
        self.x[idx, 0] = atom

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
        smiles = MolToSmiles(self.get_structure(), rootedAtAtom=root)
        smiles.replace("At", "O-")
        return "".join([(str(int(c) + ring_index) if c.isdigit() else c) for c in smiles])

    def react(self, names, types):
        """
        Override the method to the call of the reactor to modify this monomer.

        Args:
            names (List[str]): name (string representation) of the modification
            types (List[int]): Type of the parsed stings based on GlycanLexer.TYPE

        Returns:
            New monomer with the altered structure
        """
        return Reactor(self).react(names, types)

    def get_structure(self):
        """
        Compute and save the structure of this glycan.

        Returns:
            rdkit molecule representing the structure of the glycan as a graph of its non-hydrogen atoms.
        """
        if self.structure is None:
            # read the structure from the SMILES string
            self.structure = MolFromSmiles(self.smiles)

            # extract some further information from the molecule to not operate always on the molecule
            self.adjacency = GetAdjacencyMatrix(self.structure, useBO=True)
            self.ring_info = self.structure.GetRingInfo().AtomRings()
            self.x = np.zeros((self.adjacency.shape[0], 3))

            c_atoms, ringo = [], -1
            # extract some information form the molecule
            for i in range(self.adjacency.shape[0]):
                atom = self.structure.GetAtomWithIdx(i)

                # store the atom type
                self.x[i, 0] = atom.GetAtomicNum()
                if self.x[i, 0] == 6 and i in self.ring_info[0]:
                    c_atoms.append(int(i))

                # if the atom is part of any ring, store the number of that ring
                for r in range(len(self.ring_info)):
                    if i in self.ring_info[r]:
                        self.x[i, 2] = r + 1

                # identify the oxygen atom in the main ring and set its id to 10
                if self.x[i, 2] == 1 and self.x[i, 0] == 8:
                    self.x[i, 1] = 10
                    ringo = i

            self._enumerate_c_atoms(c_atoms, ringo)
        return self.structure

    def _equidistant(self, start, end):
        """
        Decider for C1 in case the previous splitting rules were all tied. This currently only fires for
        Fruf, Tagf, Sorf, Psif

        Args:
            start (int): id of the first candidate for C1 atom
            end (int): id of the second candidate for C1 atom

        Returns:
            Bool indicating that the start id is the C1 atom
        """
        c_start_candidates = np.argwhere((self.adjacency[start, :] == 1) &
                                         (self.x[:, 0] == 6) & (self.x[:, 2] == 1)).squeeze()
        c_end_candidates = np.argwhere((self.adjacency[end, :] == 1) &
                                       (self.x[:, 0] == 6) & (self.x[:, 2] == 1)).squeeze()
        if c_start_candidates.size == 1 and c_end_candidates.size == 1:
            start_ring_c = int(c_start_candidates)
            end_ring_c = int(c_end_candidates)

            start_ring_c_o_candidates = np.argwhere((self.adjacency[start_ring_c, :] == 1) &
                                                    (self.x[:, 0] == 8) & (self.x[:, 2] != 1)).squeeze()
            end_ring_c_o_candidates = np.argwhere((self.adjacency[end_ring_c, :] == 1) &
                                                  (self.x[:, 0] == 8) & (self.x[:, 2] != 1)).squeeze()

            if start_ring_c_o_candidates.size == 1 and end_ring_c_o_candidates.size == 1:
                raise UnreachableError("C1 atom cannot be detected")
            elif start_ring_c_o_candidates.size == 1:
                return True
        elif c_start_candidates.size == 1:
            return True
        return False

    def _evaluate_distance(self, start, end, ringo):
        """
        Try to decide on C1 based on their distance to the oxygen in the ring

        Args:
            start (int): id of the first candidate for C1 atom
            end (int): id of the second candidate for C1 atom
            ringo (int): index of the oxygen atom in the ring

        Returns:
            Bool indicating that the start id is the C1 atom
        """
        adj = self.adjacency.copy()

        # as we have an adjacency matrix, multiply it with itself until one of the fields is non-zero
        while adj[start, ringo] == 0 and adj[end, ringo] == 0:
            adj = adj @ self.adjacency

        # if both fields are non-zero, we cannot decide here and have to go further
        if adj[start, ringo] > 0 and adj[end, ringo] > 0:
            return self._equidistant(start, end)
        elif adj[start, ringo] > 0:
            return True
        return False

    def _enumerate_c_atoms(self, c_atoms, ringo):
        """
        Enumerate all carbon atoms starting from the first one

        Args:
            c_atoms List[int]: List of all ids of C atoms in the ring
            ringo (int): id of the oxygen atom in the ring of the monomer

        Returns:
            Nothing
        """
        # create a tree of all carbon atoms directly connected to the main ring of the monomer
        c_tree = Tree()
        stack = [(-1, c_atoms[0])]
        while len(stack) != 0:
            p_id, c_id = stack[-1]
            stack = stack[:-1]
            c_tree.add_node(c_id, p_id)

            children = np.argwhere((self.adjacency[c_id, :] == 1) & (self.x[:, 0] == 6))
            for c in children:
                if int(c) not in c_tree.nodes:
                    stack.append((c_id, int(c)))

        # find the deepest node and rehang the tree to this node
        deepest_id, _ = c_tree.deepest_node()
        c_tree = c_tree.rehang_tree(deepest_id)
        longest_c_chain = c_tree.longest_chain()

        # now the two C1 candidates can be found at the ends of the longest chain
        start, end = longest_c_chain[0], longest_c_chain[-1]

        # check conditions
        start_o_conn = np.argwhere((self.adjacency[start, :] == 1) & (self.x[:, 0] == 8) &
                                   (self.x[:, 2] != 1)).squeeze().size > 0
        end_o_conn = np.argwhere((self.adjacency[end, :] == 1) & (self.x[:, 0] == 8) &
                                 (self.x[:, 2] != 1)).squeeze().size > 0

        # decide on c1
        if start_o_conn and end_o_conn:
            if not self._evaluate_distance(start, end, ringo):
                longest_c_chain = reversed(longest_c_chain)
        elif end_o_conn:
            longest_c_chain = reversed(longest_c_chain)

        # enumerate along chain
        c_count = 0
        for c in longest_c_chain:
            c_count += 1
            self.x[c, 1] = c_count

    def find_oxygen(self, binding_c_id, check_for=None):
        """
        Find the oxygen atom that binds to the carbon atom with the provided id. The returned id may not refer to an
        oxygen atom in the ring of the monomer as this cannot bind anything. This method will report the atom id of the
        first atom type that fulfils the requirements of this method, i.e. the binding atom has to be the only one of
        its type, not in the main-ring, and connected to the given carbon atom.

        Args:
            binding_c_id (int): id of the carbon atom that participates in a binding, and we need to find the oxygen
                from
            check_for (List[int]): List if atom types as numbered in PSE to check for binding the above given carbon

        Returns:
            id referring to the atom binding the provided carbon atom and may participate in a glycan-binding.
        """
        # by default only check for binding oxygen
        if check_for is None:
            check_for = [8]

        # first find the rdkit id of the carbon atom that should bind to something
        position = np.argwhere(self.x[:, 1] == binding_c_id).squeeze()

        for check in check_for:
            # then find the candidates. There should be exactly one element in the resulting array
            candidates = np.argwhere((self.adjacency[position, :] == 1) &
                                     (self.x[:, 0] == check) & (self.x[:, 2] != 1)).squeeze()
            if candidates.size == 1:
                return int(candidates)

        raise ValueError(f"Multiple (or no) options for oxygen (or other atom type) found.")
