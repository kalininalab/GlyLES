from enum import Enum


class Monomer:
    r"""
    General interface to represent the sugar monomers in the tree.

    The resulting molecule must have some properties that are especially important for the correct implementation of
        the parser later on:
          - The C atoms have IDs 1 to 6 according to their C1-C6 naming in an glycan.
          - The O atom closing the ring has id 10, so that all other O-atoms have ID 10 higher than the C atom they are
            connected to.
          - The chirality is initial 0 and will only be set in the carbons with respect to the non-OH-group.
          - Additionally, all members of the ring of the molecule have an according flag set to true.

    Depending on the concrete implementation the actual ids of the atoms might vary, but the rest of the program will
    use the atom ids as specified here. So, an alternative naming must internally be converted to this specification.
    """

    class Config(Enum):
        r"""
        Configuration enum to represent if the monomer is in alpha, beta, or undefined configuration.
        """
        UNDEF = 0
        ALPHA = 1
        BETA = 2

    def __init__(self):
        pass

    def alpha(self):
        """
        Return this monosaccharide in its alpha conformation

        Returns:
            Monomer in alpha conformation
        """
        pass

    def beta(self):
        """
        Return this monosaccharide in its beta conformation

        Returns:
            Monomer in beta conformation
        """
        pass

    def undefined(self):
        """
        Return this monosaccharide in undefined conformation, the first carbon ring-atom will have unspecified
        chirality.

        Returns:
            Monomer in undefined conformation
        """
        pass

    def to_chirality(self, chirality):
        """
        Return this monomer in the queried chirality.

        Args:
            chirality (str): char representing the chiral conformation of the first carbon ring atom

        Returns:
            This monomer with the given (or not given) chirality at the first carbon ring atom
        """

    def get_config(self):
        """
        The current conformation relative to the first carbon ring-atom, i.e. alpha, beta or unspecified.

        Returns:
            Config-Tag according to the conformation this monomer represents
        """
        pass

    def is_non_chiral(self):
        """
        Check if this monomer represents a non-chiral molecule

        Returns:
            boolean indicating chirality of this monomer
        """
        pass

    def get_dummy_atoms(self):
        """
        Specify some dummy atoms that are used to mark oxygen atoms that will participate in bindings between glycans.

        Returns:
            Two lists:
                * one with atoms that are given as the "atom" argument to the mark-method to replace the oxygen atoms
                * the string representation of the atoms from above, i.e. how the atoms above will be represented in a
                  SMILES string
        """
        pass

    def root_atom_id(self, binding_c_id):
        """
        Get ID of atom atom that will bind the parent monomer in the glycan. This ID will be given as root argument to
        the to_smiles method.

        Args:
            binding_c_id (int): Integer at which c-position this monomer binds its parent

        Returns:
            id of the atom that binds to the parent
        """
        pass

    def mark(self, position, atom):
        """
        Mark the oxygen atom linked to the carbon atom at the given position ready to participate in the bounding.

        Args:
            position: id of the carbon atom whose oxygen atom will from the binding
            atom: atom to replace the binding oxygen with

        Returns:
            Nothing
        """
        pass

    def to_smiles(self, root, ring_index):
        """
        Convert the this monomer into a SMILES string representation.

        Args:
            root (int): index of the root atom
            ring_index: index of the rings in the atom

        Returns:
            SMILES string representation of this molecule
        """
        pass

    @staticmethod
    def from_string(mono):
        """
        Get a monomer from the string representation of the monomer.

        Args:
            mono: name of the monomer to get the object from

        Returns:
            An instance of this class to be used to build up the complete polysaccharide
        """
        pass
