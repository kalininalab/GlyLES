import numpy as np
from rdkit.Chem import MolFromSmiles, MolToSmiles, GetAdjacencyMatrix

from glyles.glycans.enum_c import enumerate_carbon
from glyles.glycans.reactor import SMILESReaktor
from glyles.glycans.utils import Config, find_isomorphism_nx
from glyles.grammar.GlycanLexer import GlycanLexer


class Monomer:

    def __init__(self, origin=None, **kwargs):
        """
        Initialize the monomer using the super method. Additionally, some fields are initialized to describe the
        structure of the monomer according to the specification of the monomer-parent class

        Args:
            origin (Monomer): Other monomer to use to initialize this object
            **kwargs: arguments to initialize monomer if object is None.
        """
        if isinstance(origin, tuple):
            origin = origin[0]

        if origin is None:
            self.name = kwargs["name"]
            self.smiles = kwargs["smiles"]
            self.__root_smiles = kwargs["smiles"]
            self.structure = kwargs.get("struct", None)
            self.config = kwargs["config"]
            self.isomer = kwargs["isomer"]
            self.lactole = kwargs["lactole"]
            self.recipe = kwargs["recipe"]
            self.c1_find = kwargs.get("c1_find", None)
        else:
            self.name = origin.get_name()
            self.smiles = origin.get_smiles()
            self.__root_smiles = origin.root_smiles()
            self.c1_find = origin.get_c1_finder()
            self.structure = origin.get_structure()
            self.config = origin.get_config()
            self.isomer = origin.get_isomer()
            self.lactole = origin.get_lactole()
            self.recipe = origin.get_recipe()
            self.adjacency = origin.get_adjacency()
            self.ring_info = origin.get_ring_info()
            self.x = origin.get_features()

    def root_smiles(self):
        """
        Access the SMILES string of the root monomer. This is before any functional group is attached.
        It is as described in the factories.

        Returns:
            SMILES string of the root monomer of this monosaccharide
        """
        return self.__root_smiles

    def get_name(self, full=False):
        """
        Returns the name of this monomer as three-letter code (eventually longer for more fancy monosaccharides with
        more complex side chains).

        Args:
            full (bool): Flag indicating to also output recipe of monomer

        Returns:
            The name of this monomer
        """
        if full:
            return "".join((s for s, _ in self.recipe))
        return self.name

    def get_smiles(self):
        """
        Returns the SMILES representation of this monomer. Attention: These methods return the SMILES that is used to
        initialize this monomer. This is different from the to_smiles method of this class that returns the SMILES
        string with added place-holders that is used for the generation of the SMILES representation of the complete
        glycan

        Returns:
            The SMILES string that was used for initialization of this monomer
        """
        return self.smiles

    def get_c1_finder(self):
        """
        This returns a method used to identify the C1 atom in an additionally provided molecule object.

        Returns:
            The callable object that, given the structure, finds the RDKit ID of C1
        """
        return self.c1_find

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
        return Monomer(factory.create(recipe))

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
        return Monomer(factory.create(recipe))

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
        return Monomer(factory.create(recipe))

    def to_chirality(self, chirality, factory):
        """
        Return this monomer in the queried chirality.

        Args:
            chirality (str): char representing the chiral conformation of the first carbon ring atom
            factory (MonomerFactory): factory instance to be used to generate the monomers

        Returns:
            This monomer with the given (or not given) chirality at the first carbon ring atom
        """
        chirality = chirality.lower()
        if chirality == "a":
            return self.alpha(factory)
        if chirality == "b":
            return self.beta(factory)
        return self.undefined(factory)

    def get_config(self):
        """
        The current conformation relative to the first carbon ring-atom, i.e. alpha, beta or unspecified.

        Returns:
            Config-Tag according to the conformation this monomer represents
        """
        return self.config

    def get_isomer(self):
        """
        The current enantiomer of this monomer, i.e. L-form or D-form.

        Returns:
            Enantiomer-Tag according to the isomer this monomer represents
        """
        return self.isomer

    def get_lactole(self):
        """
        The current lactole-form of this monomer, i.e. if it's a 5-ring or a 6-ring molecule.

        Returns:
            Lactole-Tag according to the lactole-form this monomer represents
        """
        return self.lactole

    def get_recipe(self):
        """
        Return the recipe on how this monomer was created

        Returns:
            Recipe as list of modifications, confirmations and the root monomer as type List[Tuple[str, int]]
        """
        return self.recipe

    def is_non_chiral(self):
        """
        Check if this monomer represents a non-chiral molecule

        Returns:
            boolean indicating chirality of this monomer
        """
        return self.config == Config.UNDEF

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

    @staticmethod
    def get_dummy_atoms():
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
        return [34, 52, 84, 85], ["[SeH]", "[TeH]", "[PoH]", "[At]"]

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

    def to_smiles(self, ring_index, root_idx=None, root_id=None):
        """
        Convert this monomer into a SMILES string representation.
        Use the implementation of the SMILES algorithm fitted to the needs of glycans.

        Args:
            ring_index (int): index of the rings in the atom
            root_idx (int): index of the root atom
            root_id (int): RDKit ID of root atom

        Returns:
            SMILES string representation of this molecule
        """
        assert root_idx is not None or root_id is not None, "Either Index or ID has to be provided"
        if root_id is None:
            if np.where(self.x[:, 1] == root_idx)[0].size != 0:
                root_id = int(np.where(self.x[:, 1] == root_idx)[0])
            else:
                root_id = int(np.where(self.x[:, 1] == 1)[0])
        smiles = MolToSmiles(self.get_structure(), rootedAtAtom=root_id)
        return "".join([((f"%{int(c) + ring_index}" if int(c) + ring_index >= 10
                          else f"{int(c) + ring_index}") if c.isdigit() else c) for c in smiles])

    def react(self, names, types):
        """
        Override the method to the call of the reactor to modify this monomer.

        Args:
            names (List[str]): name (string representation) of the modification
            types (List[int]): Type of the parsed stings based on GlycanLexer.TYPE

        Returns:
            New monomer with the altered structure
        """
        return SMILESReaktor(self).react(names, types)

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
            rings = self.structure.GetRingInfo().AtomRings()
            if len(rings) > 0 and self.name != "Inositol":
                self.ring_info = [None]
                for ring in rings:
                    found_ox = False
                    for atom in ring:
                        if self.structure.GetAtomWithIdx(atom).GetAtomicNum() == 8:
                            self.ring_info[0] = ring
                            found_ox = True
                            break
                    if not found_ox:
                        self.ring_info.append(ring)
            else:
                self.ring_info = rings
            self.x = np.zeros((self.adjacency.shape[0], 4), dtype=int)

            # extract some information form the molecule
            for i in range(self.adjacency.shape[0]):
                atom = self.structure.GetAtomWithIdx(i)

                # store the atom type
                self.x[i, 0] = atom.GetAtomicNum()

                # if the atom is part of any ring, store the number of that ring
                for r in range(len(self.ring_info)):
                    if self.ring_info[r] is not None and i in self.ring_info[r]:
                        self.x[i, 2] = r + 1

                # identify the oxygen atom in the main ring and set its id to 100
                if self.x[i, 2] == 1 and self.x[i, 0] == 8:
                    self.x[i, 1] = 100

            # identify isomorphic atoms. The reference for the isomorphism test is the root SMILES
            iso = list(find_isomorphism_nx(self.smiles, self.root_smiles(), self.name, self.c1_find).keys())
            if len(iso) == 0:
                iso = list(find_isomorphism_nx(self.root_smiles(), self.smiles, self.name, self.c1_find).values())
            if len(iso) == 0:
                self.x[:, 3] = 1
            else:
                self.x[iso, 3] = 1

            # Enumerate all carbons
            enumerate_carbon(self)

        return self.structure

    def find_oxygen(self, binding_c_id):
        """
        Find the oxygen atom that binds to the carbon atom with the provided id. The returned id may not refer to an
        oxygen atom in the ring of the monomer as this cannot bind anything. This method will report the atom id of the
        first atom type that fulfils the requirements of this method, i.e. the binding atom has to be the only one of
        its type, not in the main-ring, and connected to the given carbon atom.

        Args:
            binding_c_id (int): id of the carbon atom that participates in a binding, and we need to find the oxygen
                from

        Returns:
            The RDKit-ID referring to the atom binding the provided carbon atom and may participate in a glycan-binding.
            In case the carbon is bound to neither an oxygen nor a nitrogen, the RDKit-ID of the carbon is returned.
        """
        # first find the rdkit id of the carbon atom that should bind to something
        position = np.argwhere(self.x[:, 1] == binding_c_id).squeeze()

        multiple = False
        for check in [8, 7]:
            # then find the candidates. There should be exactly one element in the resulting array
            candidates = np.argwhere((self.adjacency[position, :] == 1) &
                                     (self.x[:, 0] == check) & (self.x[:, 2] != 1)).squeeze()
            if candidates.size == 1:
                return int(candidates)
            elif candidates.size > 0:
                multiple = True

        if not multiple and position.size == 1:
            return int(position)

        raise ValueError(f"Multiple options for oxygen (or other atom type) found.")
