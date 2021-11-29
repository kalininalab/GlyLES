import numpy as np
from rdkit.Chem import MolFromSmiles, MolToSmiles, GetAdjacencyMatrix

from glyles.glycans.monomer import Monomer


# TODO: Extract structure generation and annotation into script and save a file for each structure
# reduce time for preprocessing


class RDKitMonomer(Monomer):
    def __init__(self, origin=None, **kwargs):
        """
        Initialize the monomer using the super method. Additionally, some fields are initialized to describe the
        structure of the monomer according to the specification of the monomer-parent class

        Args:
            origin (Monomer): Other monomer to use to initialize this object
            **kwargs: arguments to initialize monomer if object is None. Must include name, smiles, and config
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
        Get ID of atom atom that will bind the parent monomer in the glycan. This ID will be given as root argument to
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
        Convert the this monomer into a SMILES string representation.
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
            # read the structure from the smiles string
            self._structure = MolFromSmiles(self._smiles)

            # extract some further information from the molecule to not operate always on the molecule
            self._adjacency = GetAdjacencyMatrix(self._structure)
            self._ring_info = self._structure.GetRingInfo().AtomRings()
            self._x = np.zeros((self._adjacency.shape[0], 3))

            # extract some information form the molecule
            for i in range(self._adjacency.shape[0]):
                atom = self._structure.GetAtomWithIdx(i)

                # store the atom type
                self._x[i, 0] = atom.GetAtomicNum()

                # if the atom is part of any ring, store the number of that ring
                for r in range(len(self._ring_info)):
                    if i in self._ring_info[r]:
                        self._x[i, 2] = r + 1

                # identify the oxygen atom in the main ring and set its id to 10
                if self._x[i, 2] == 1 and self._x[i, 0] == 8:
                    self._x[i, 1] = 10

            # check in which order the ring has to be trversed to find C1 as first carbon atom behind the ring-oxygen
            main_ring = list(self._ring_info[0]) + list(self._ring_info[0])
            if not self.__clockwise():
                main_ring = list(reversed(main_ring))

            # TODO: This only works if the C1 atom is part of the ring! (i.e. not for fructose)
            # for all carbon atoms in the main ring, set its id according to the distance to the oxygen atom
            # iterate clockwise, i.e. from ring-O to C1/C2 and so on ...
            carbon_count = 1
            # Insert code for c1-atoms outside of the ring here
            o_index = main_ring.index(np.argwhere((self._x[:, 0] == 8) & (self._x[:, 2] == 1)))
            for r in self._ring_info[0]:
                if self._x[r, 0] != 8:
                    self._x[r, 1] = carbon_count + main_ring.index(r, o_index) - o_index - 1

            # assign ids to the remaining carbon atoms
            highest_index = np.max(self._x[:, 1][np.argwhere(self._x[:, 0] == 6)])
            # c6_candidates = np.argwhere((self._x[:, 1] == 6) & (self._x[:, 2] == 0))

            c5 = np.argwhere((self._x[:, 1] == 5).squeeze())
            c6 = np.argwhere(((self._x[:, 0] == 6) & (self._x[:, 2] == 0) & (self._adjacency[c5, :] == 1)).squeeze())
            self._x[c6, 1] = highest_index + 1
        return self._structure

    def __clockwise(self):
        """
        Check if the main ring in the monomer is iterated in clockwise manner. RDKit assigns the atom_ids in the order
        it is read in during parsing the monomer SMILES strings.

        Returns:
            True if the main ring is traversed clockwise when iterating over all ring member in order of their RDKit-ID
        """
        position = np.argwhere(self._x[:, 1] == 10)

        # then find the candidates. There should be exactly one element in the resulting array
        candidates = np.argwhere((self._adjacency[position, :] == 1).squeeze())
        if len(candidates) != 2:
            raise ValueError("Parsing exception")

        # check if the candidate with the lower id fulfills the criteria to be the C1 atom
        first_c_cons = np.argwhere((self._adjacency[candidates[0], :] == 1).squeeze())
        if np.sum(self._x[first_c_cons, 0] == 6) == 1:
            o_index = self._ring_info[0].index(np.argwhere((self._x[:, 0] == 8) & (self._x[:, 2] == 1)))
            c_index = self._ring_info[0].index(candidates[0])
            if o_index < c_index:
                return True
        return False

    def __find_oxygen(self, binding_c_id):
        """
        Find the oxygen atom that binds to the carbon atom with the provided id. The returned id may not refer to an
        oxygen atom in the ring of the monomer as this cannot bind anything.

        Args:
            binding_c_id (int): id of the carbon atom that participates in a binding and we need the oxygen from

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
