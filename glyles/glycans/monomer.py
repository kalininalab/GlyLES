# from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.utils import Config, Lactole, Enantiomer


class Monomer:
    """
    General interface to represent the sugar monomers in the tree.

    The resulting molecule must have some properties that are especially important for the correct implementation of
        the parser later on:
          - The C atoms have IDs 1 to 6 according to their C1-C6 naming in a glycan.
          - The O atom closing the ring has id 10, so that all other O-atoms have ID 10 higher than the C atom they are
            connected to.
          - The chirality is initial 0 and will only be set in the carbons with respect to the non-OH-group.
          - Additionally, all members of the ring of the molecule have an according flag set to true.

    Depending on the concrete implementation the actual ids of the atoms might vary, but the rest of the program will
    use the atom ids as specified here. So, an alternative naming must internally be converted to this specification.
    """

    def __init__(self, origin=None, **kwargs):
        """
        Initialize the Monomer with the provided arguments or from another monomer

        Args:
            origin (Monomer): Other monomer to use to initialize this object
            **kwargs: arguments to initialize monomer if object is None. Must include name, SMILES, and config
        """
        if origin is None:
            self._name = kwargs["name"]
            self._smiles = kwargs["smiles"]
            self._structure = kwargs.get("struct", None)
            self._config = kwargs["config"]
            self._isomer = kwargs["isomer"]
            self._lactole = kwargs["lactole"]
        else:
            self._name = origin.get_name()
            self._smiles = origin.get_smiles()
            self._structure = origin.get_structure()
            self._config = origin.get_config()
            self._isomer = origin.get_isomer()
            self._lactole = origin.get_lactole()

    def get_name(self):
        """
        Returns the name of this monomer as three-letter code (eventually longer for more fancy monosaccharides with
        more complex side chains).

        Returns:
            The name of this monomer
        """
        return self._name

    def get_smiles(self):
        """
        Returns the SMILES representation of this monomer. Attention: These methods return the SMILES that is used to
        initialize this monomer. This is different from the to_smiles method of this class that returns the SMILES
        string with added place-holders that is used for the generation of the SMILES representation of the complete
        glycan

        Returns:
            The SMILES string that was used for initialization of this monomer
        """
        return self._smiles

    def get_structure(self):
        """
        Returns the structure of this monomer. This might be Null in case the structure wasn't inferred from the
        SMILES so far

        Returns:
            Structure depending on the methods chosen to represent the structure of this monomer, might be Null
        """
        return self._structure

    def alpha(self, factory):
        """
        Return this monosaccharide in its alpha conformation.

        Args:
            factory (MonomerFactory): factory instance to be used to generate the monomers

        Returns:
            Monomer in alpha conformation
        """
        return Monomer(factory["A_" + self._name])

    def beta(self, factory):
        """
        Return this monosaccharide in its beta conformation.

        Args:
            factory (MonomerFactory): factory instance to be used to generate the monomers

        Returns:
            Monomer in beta conformation
        """
        return Monomer(factory["B_" + self._name])

    def undefined(self, factory):
        """
        Return this monosaccharide in undefined conformation, the first carbon ring-atom will have unspecified.
        chirality.

        Args:
            factory (MonomerFactory): factory instance to be used to generate the monomers

        Returns:
            Monomer in undefined conformation
        """
        return Monomer(factory[self._name])

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
        return self._config

    def get_isomer(self):
        """
        The current enantiomer of this monomer, i.e. L-form or D-form.

        Returns:
            Enantiomer-Tag according to the isomer this monomer represents
        """
        return self._isomer

    def get_lactole(self):
        """
        The current lactole-form of this monomer, i.e. if it's a 5-ring or a 6-ring molecule.

        Returns:
            Lactole-Tag according to the lactole-form this monomer represents
        """
        return self._lactole

    def is_non_chiral(self):
        """
        Check if this monomer represents a non-chiral molecule

        Returns:
            boolean indicating chirality of this monomer
        """
        return self._config == Config.UNDEF

    def get_dummy_atoms(self):
        """
        Specify some dummy atoms that are used to mark oxygen atoms that will participate in bindings between glycans.

        Returns:
            Two lists:
                * one with atoms that are given as the "atom" argument to the mark-method to replace the oxygen atoms
                * the string representation of the atoms from above, i.e. how the atoms above will be represented in a
                  SMILES string
        """
        raise NotImplementedError

    def root_atom_id(self, binding_c_id):
        """
        Get ID of atom that will bind the parent monomer in the glycan. This ID will be given as root argument to
        the to_smiles method.

        Args:
            binding_c_id (int): Integer at which c-position this monomer binds its parent

        Returns:
            id of the atom that binds to the parent
        """
        raise NotImplementedError

    def mark(self, position, atom):
        """
        Mark the oxygen atom linked to the carbon atom at the given position ready to participate in the bounding.

        Args:
            position (int): id of the carbon atom whose oxygen atom will from the binding
            atom (object): atom to replace the binding oxygen with

        Returns:
            Nothing
        """
        raise NotImplementedError

    def to_smiles(self, root, ring_index):
        """
        Convert this monomer into a SMILES string representation.

        Args:
            root (int): index of the root atom
            ring_index (int): index of the rings in the atom

        Returns:
            SMILES string representation of this molecule
        """
        raise NotImplementedError

    def __get_structure(self):
        """
        Compute and save the structure of this glycan.

        Returns:
            Representation the structure of the glycan as a graph.
        """
        raise NotImplementedError
