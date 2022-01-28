import warnings

import numpy as np
from rdkit.Chem import MolFromSmiles, MolToSmiles, GetAdjacencyMatrix
from rdkit.Chem.rdchem import Atom, EditableMol, BondType

from glyles.glycans.monomer import Monomer
from glyles.glycans.utils import UnreachableError, Tree, ModificationNotImplementedWarning
from glyles.grammar.GlycanLexer import GlycanLexer


def not_implemented_warning(mod):
    print(f"ModificationNotImplementedWarning: {mod} Modification not implemented. The returned molecule will not have "
          f"this modification")
    warnings.warn(f"{mod} Modification not implemented. The returned molecule will not have "
                  f"this modification", ModificationNotImplementedWarning)


class RDKitMonomer(Monomer):
    class Reactor:
        """
        Class with access to protected classes fields managing the modifications of a root monomer.
        """

        def __init__(self, monomer):
            """
            Initialize the reactor with the monomer to modify.

            Args:
                monomer (RDKitMonomer): monomer to be modified
            """
            self.monomer = monomer
            if self.monomer._structure is None:
                self.monomer._get_structure()

        def react(self, names, types):
            """
            Manage the parsed modifications and apply them in turn.

            Args:
                names (List[str]): name (string representation) of the modification
                types (List[int]): Type of the parsed stings based on GlycanLexer.TYPE

            Returns:
                Modified monomer
            """
            for name, t in zip(names, types):
                if t != GlycanLexer.MOD:
                    continue
                if len(name) == 1:
                    if name[0] == "N":
                        self.set_nitrogen()
                    if name[0] == "A":
                        self.make_acid()
                elif len(name) == 2:
                    if name[0].isdigit():
                        if name[1] == "d":  # ?d
                            # raise NotImplementedError("Deoxygenation not implemented yet")
                            not_implemented_warning(name)
                        elif name[1] == "e":  # ?d
                            # raise NotImplementedError("Deoxygenation not implemented yet")
                            not_implemented_warning(name)
                        elif name[1] == "S":
                            self.add_sulfur(int(name[0]))
                        elif name[1] == "P":
                            self.add_phosphate(int(name[0]))
                        elif name[1] == "N":
                            self.set_nitrogen(position=int(name[0]))
                    else:
                        if name == "Ac":
                            self.add_acid(position=5)
                elif len(name) == 3:
                    if name == "NAc":
                        self.add_acid(pos=self.set_nitrogen())
                    elif name[0].isdigit():
                        if name.endswith("Me"):
                            self.add_methyl(position=int(name[0]))
                        elif name.endswith("Ac"):
                            self.add_acid(position=int(name[0]))
                elif len(name) == 7:
                    if name.endswith("Me-"):
                        self.add_methyl(position=int(name[0]))
                    if name.endswith("Ac-"):
                        self.add_acid(position=int(name[0]))

            return self.monomer

        def add_sulfur(self, position):
            """
            Add a SO3- group at the oxygen of the specified position.
            Example: Gal -> Gal3S

            Args:
                position (int): id of the carbon where to add the SO3- group to the bound oxygen

            Returns:
                Nothing
            """
            pos = self.monomer._find_oxygen(position)
            emol = EditableMol(self.monomer._structure)

            s_id = EditableMol.AddAtom(emol, Atom(16))
            o1_id = EditableMol.AddAtom(emol, Atom(8))
            o2_id = EditableMol.AddAtom(emol, Atom(8))
            o3_id = EditableMol.AddAtom(emol, Atom(8))
            EditableMol.AddBond(emol, s_id, o1_id, order=BondType.DOUBLE)
            EditableMol.AddBond(emol, s_id, o2_id, order=BondType.SINGLE)
            EditableMol.AddBond(emol, s_id, o3_id, order=BondType.DOUBLE)
            EditableMol.AddBond(emol, pos, s_id, order=BondType.SINGLE)

            self.monomer._structure = emol.GetMol()

            new_x, new_adj = self._extend_matrices(4)
            new_x[s_id:, 0] = [16, 8, 8, 8]
            self.monomer._x = new_x

            new_adj[s_id, o1_id] = 1
            new_adj[s_id, o2_id] = 1
            new_adj[s_id, o3_id] = 1
            new_adj[s_id, pos] = 1
            new_adj[o1_id, s_id] = 1
            new_adj[o2_id, s_id] = 1
            new_adj[o3_id, s_id] = 1
            new_adj[pos, s_id] = 1

            self.monomer._adjacency = new_adj

        def add_phosphate(self, position):
            """
            Add a PO3 group at the oxygen of the specified position.
            Example: Gal -> Gal3P

            Args:
                position (int): id of the carbon where to add the PO3 group to the bound oxygen

            Returns:
                Nothing
            """
            pos = self.monomer._find_oxygen(position)
            emol = EditableMol(self.monomer._structure)

            s_id = EditableMol.AddAtom(emol, Atom(15))
            o1_id = EditableMol.AddAtom(emol, Atom(8))
            o2_id = EditableMol.AddAtom(emol, Atom(8))
            o3_id = EditableMol.AddAtom(emol, Atom(8))
            EditableMol.AddBond(emol, s_id, o1_id, order=BondType.SINGLE)
            EditableMol.AddBond(emol, s_id, o2_id, order=BondType.SINGLE)
            EditableMol.AddBond(emol, s_id, o3_id, order=BondType.DOUBLE)
            EditableMol.AddBond(emol, pos, s_id, order=BondType.SINGLE)

            self.monomer._structure = emol.GetMol()

            new_x, new_adj = self._extend_matrices(4)
            new_x[s_id:, 0] = [15, 8, 8, 8]
            self.monomer._x = new_x

            new_adj[s_id, o1_id] = 1
            new_adj[s_id, o2_id] = 1
            new_adj[s_id, o3_id] = 1
            new_adj[s_id, pos] = 1
            new_adj[o1_id, s_id] = 1
            new_adj[o2_id, s_id] = 1
            new_adj[o3_id, s_id] = 1
            new_adj[pos, s_id] = 1

            self.monomer._adjacency = new_adj

        def add_acid(self, position=None, pos=None):
            """
            Add an acid group to a specific position. Here the position can be provided either as the C-index
            (position) or as the rdkit id of the atom where to append the acid group (implemented for NAc). Exactly one
            of both arguments must be provided.
            Example: GalN -> GalNAc or Gal -> Gal5Ac

            Args:
                position (int): index of the c-atom where to append the acid group
                pos (int): rdkit id of the atom where to append the acid group

            Returns:
                Nothing
            """
            if (position is None) == (pos is None):
                raise ValueError()

            if position is not None:
                pos = self.monomer._find_oxygen(position, check_for=[8, 7])

            emol = EditableMol(self.monomer._structure)

            c1_id = EditableMol.AddAtom(emol, Atom(6))
            c2_id = EditableMol.AddAtom(emol, Atom(6))
            o1_id = EditableMol.AddAtom(emol, Atom(8))
            EditableMol.AddBond(emol, c1_id, pos, order=BondType.SINGLE)
            EditableMol.AddBond(emol, c2_id, c1_id, order=BondType.SINGLE)
            EditableMol.AddBond(emol, o1_id, c1_id, order=BondType.DOUBLE)

            self.monomer._structure = emol.GetMol()

            new_x, new_adj = self._extend_matrices(3)
            new_x[c1_id:, 0] = [6, 6, 8]
            self.monomer._x = new_x

            new_adj[c1_id, c2_id] = 1
            new_adj[c1_id, o1_id] = 1
            new_adj[c1_id, pos] = 1
            new_adj[c2_id, c1_id] = 1
            new_adj[o1_id, c1_id] = 1
            new_adj[pos, c1_id] = 1

            self.monomer._adjacency = new_adj

        def add_methyl(self, position):
            """

            Args:
                position:

            Returns:

            """
            pos = self.monomer._find_oxygen(position)

            emol = EditableMol(self.monomer._structure)

            c_id = EditableMol.AddAtom(emol, Atom(6))
            EditableMol.AddBond(emol, c_id, pos, order=BondType.SINGLE)

            self.monomer._structure = emol.GetMol()

            new_x, new_adj = self._extend_matrices(1)
            new_x[c_id, 0] = 6
            self.monomer._x = new_x

            new_adj[c_id, pos] = 1
            new_adj[pos, c_id] = 1

            self.monomer._adjacency = new_adj

        def set_nitrogen(self, position=2):
            """
            Change an oxygen to a nitrogen.
            Example: Gal -> GalN

            Args:
                position (int): position of a carbon atom at which the bound o atom should be replaced by N

            Returns:
                rdkit id of the atom that is now a nitrogen
            """
            pos = self.monomer._find_oxygen(position, check_for=[8, 7])
            self.monomer._structure.GetAtomWithIdx(pos).SetAtomicNum(7)
            self.monomer._x[pos, 0] = 7
            return pos

        def make_acid(self):
            """

            Returns:

            """
            emol = EditableMol(self.monomer._structure)

            c_id = int(np.argwhere(self.monomer._x[:, 1] == 6))
            o_id = EditableMol.AddAtom(emol, Atom(8))
            EditableMol.AddBond(emol, o_id, c_id, order=BondType.DOUBLE)

            self.monomer._structure = emol.GetMol()

            new_x, new_adj = self._extend_matrices(1)
            new_x[o_id, 0] = 8
            self.monomer._x = new_x

            new_adj[c_id, o_id] = 1
            new_adj[o_id, c_id] = 1
            self.monomer._adjacency = new_adj

        def _extend_matrices(self, count):
            """

            Args:
                count:

            Returns:

            """
            tmp_x = np.zeros((self.monomer._x.shape[0] + count, self.monomer._x.shape[1]))
            tmp_x[:self.monomer._x.shape[0], :] = self.monomer._x
            tmp_adj = np.zeros((self.monomer._adjacency.shape[0] + count, self.monomer._adjacency.shape[1] + count))
            tmp_adj[:self.monomer._adjacency.shape[0], :self.monomer._adjacency.shape[1]] = self.monomer._adjacency
            return tmp_x, tmp_adj

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
            self._get_structure()

    def alpha(self, factory):
        """
        Return this monosaccharide in its alpha conformation.

        Args:
            factory (MonomerFactory): factory instance to be used to generate the monomers

        Returns:
            Monomer in alpha conformation
        """
        recipe = [(v, t) for v, t in self._recipe if t != GlycanLexer.TYPE]
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
        recipe = [(v, t) for v, t in self._recipe if t != GlycanLexer.TYPE]
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
        recipe = [(v, t) for v, t in self._recipe if t != GlycanLexer.TYPE]
        return RDKitMonomer(factory.create(recipe))

    def get_adjacency(self):
        """
        Get the adjacency matrix of the atoms in this monomer.

        Returns:
            Adjacency matrix of all non-hydrogen atoms in this monomer
        """
        if self._adjacency is None:
            self._get_structure()
        return self._adjacency

    def get_ring_info(self):
        """
        Get information of all atom-ids in the rings of this monomer.

        Returns:
            Tuple of tuples with the atom-ids from rdkit in the monomer
        """
        if self._ring_info is None:
            self._get_structure()
        return self._ring_info

    def get_features(self):
        """
        Get a feature matrix from this monomer. The features contain information about atom type, ids according to the
        specification in monomer.py, and ring memberships.

        Returns:
            A numpy array of shape Nx3 containing the extracted features for all atoms in this molecule
        """
        if self._x is None:
            self._get_structure()
        return self._x

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
        return self._find_oxygen(binding_c_id)

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
        idx = self._find_oxygen(position)
        self._get_structure().GetAtomWithIdx(idx).SetAtomicNum(atom)
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
        smiles = MolToSmiles(self._get_structure(), rootedAtAtom=root)
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
        return RDKitMonomer.Reactor(self).react(names, types)

    def _get_structure(self):
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

            self._enumerate_c_atoms(c_atoms, ringo)
        return self._structure

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
        c_start_candidates = np.argwhere((self._adjacency[start, :] == 1) &
                                         (self._x[:, 0] == 6) & (self._x[:, 2] == 1)).squeeze()
        c_end_candidates = np.argwhere((self._adjacency[end, :] == 1) &
                                       (self._x[:, 0] == 6) & (self._x[:, 2] == 1)).squeeze()
        if c_start_candidates.size == 1 and c_end_candidates.size == 1:
            start_ring_c = int(c_start_candidates)
            end_ring_c = int(c_end_candidates)

            start_ring_c_o_candidates = np.argwhere((self._adjacency[start_ring_c, :] == 1) &
                                                    (self._x[:, 0] == 8) & (self._x[:, 2] != 1)).squeeze()
            end_ring_c_o_candidates = np.argwhere((self._adjacency[end_ring_c, :] == 1) &
                                                  (self._x[:, 0] == 8) & (self._x[:, 2] != 1)).squeeze()

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
        adj = self._adjacency.copy()

        # as we have an adjacency matrix, multiply it with itself until one of the fields is non-zero
        while adj[start, ringo] == 0 and adj[end, ringo] == 0:
            adj = adj @ self._adjacency

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

            children = np.argwhere(self._adjacency[c_id] & (self._x[:, 0] == 6))
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
        start_o_conn = np.argwhere((self._adjacency[start, :] == 1) & (self._x[:, 0] == 8) &
                                   (self._x[:, 2] != 1)).squeeze().size > 0
        end_o_conn = np.argwhere((self._adjacency[end, :] == 1) & (self._x[:, 0] == 8) &
                                 (self._x[:, 2] != 1)).squeeze().size > 0

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
            self._x[c, 1] = c_count

    def _find_oxygen(self, binding_c_id, check_for=None):
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
        position = np.argwhere(self._x[:, 1] == binding_c_id).squeeze()

        for check in check_for:
            # then find the candidates. There should be exactly one element in the resulting array
            candidates = np.argwhere((self._adjacency[position, :] == 1) &
                                     (self._x[:, 0] == check) & (self._x[:, 2] != 1)).squeeze()
            if candidates.size == 1:
                return int(candidates)

        raise ValueError(f"Multiple (or no) options for oxygen (or other atom type) found.")
