from enum import Enum


class Monomer:
    """
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
        """
        Configuration enum to represent if the monomer is in alpha, beta, or undefined configuration.
        """
        UNDEF = 0
        ALPHA = 1
        BETA = 2

    def __init__(self, origin=None, **kwargs):
        """
        Initialize the Monomer with the provided arguments or from another monomer

        Args:
            origin (Monomer): Other monomer to use to initialize this object
            **kwargs: arguments to initialize monomer if object is None. Must include name, smiles, and config
        """
        if origin is None:
            self._name = kwargs["name"]
            self._smiles = kwargs["smiles"]
            self._structure = kwargs.get("struct", None)
            self._config = kwargs.get("config", Monomer.Config.UNDEF)
        else:
            self._name = origin.get_name()
            self._smiles = origin.get_smiles()
            self._structure = origin.get_structure()
            self._config = origin.get_config()

    def get_name(self):
        """
        Returns the name of this monomer as three letter code (eventually longer for more fancy monosaccharides with
        more complex side chains).

        Returns:
            The name of this monomer
        """
        return self._name

    def get_smiles(self):
        """
        Returns the smiles representation of this monomer. Attention: This methods returns the SMILES that is used to
        initialize this monomer. This is different from the to_smiles method of this class that returns the SMILES
        string with added place-holders that is used for the generation of the smiles representation of the complete
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

    def alpha(self):
        """
        Return this monosaccharide in its alpha conformation.

        Returns:
            Monomer in alpha conformation
        """
        return self.from_string("A" + self._name[-3:])

    def beta(self):
        """
        Return this monosaccharide in its beta conformation.

        Returns:
            Monomer in beta conformation
        """
        return self.from_string("B" + self._name[-3:])

    def undefined(self):
        """
        Return this monosaccharide in undefined conformation, the first carbon ring-atom will have unspecified.
        chirality.

        Returns:
            Monomer in undefined conformation
        """
        return self.from_string(self._name[-3:])

    def to_chirality(self, chirality):
        """
        Return this monomer in the queried chirality.
        Args:
            chirality (str): char representing the chiral conformation of the first carbon ring atom

        Returns:
            This monomer with the given (or not given) chirality at the first carbon ring atom
        """
        chirality = chirality.lower()
        if chirality == "a":
            return self.alpha()
        if chirality == "b":
            return self.beta()
        return self.undefined()

    def get_config(self):
        """
        The current conformation relative to the first carbon ring-atom, i.e. alpha, beta or unspecified.

        Returns:
            Config-Tag according to the conformation this monomer represents
        """
        return self._config

    def is_non_chiral(self):
        """
        Check if this monomer represents a non-chiral molecule

        Returns:
            boolean indicating chirality of this monomer
        """
        return self._config == Monomer.Config.UNDEF

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
        Get ID of atom atom that will bind the parent monomer in the glycan. This ID will be given as root argument to
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
        Convert the this monomer into a SMILES string representation.

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

    @staticmethod
    def from_string(mono):
        """
        Get an instance of a glycan according to the provided string representation of the glycan.

        Args:
            mono (str): string representation of the glycan of interest

        Returns:
            Glycan according to the monosaccharide provided via mono
        """
        return {
            "FUC": Monomer(name="Fuc", config=Monomer.Config.UNDEF,
                           smiles="C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O"),
            "AFUC": Monomer(name="Fuc", config=Monomer.Config.ALPHA,
                            smiles="C[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O",),
            "BFUC": Monomer(name="Fuc", config=Monomer.Config.BETA,
                            smiles="C[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"),

            "GAL": Monomer(name="Gal", config=Monomer.Config.UNDEF,
                           smiles="OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"),
            "AGAL": Monomer(name="Gal", config=Monomer.Config.ALPHA,
                            smiles="OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"),
            "BGAL": Monomer(name="Gal", config=Monomer.Config.BETA,
                            smiles="OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"),

            "GAL3S": Monomer(name="Gal6S", config=Monomer.Config.UNDEF,
                             smiles="O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O"),
            "AGAL3S": Monomer(name="Gal6S", config=Monomer.Config.ALPHA,
                              smiles="O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1O"),
            "BGAL3S": Monomer(name="Gal6S", config=Monomer.Config.BETA,
                              smiles="O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O"),

            "GAL6S": Monomer(name="Gal6S", config=Monomer.Config.UNDEF,
                             smiles="O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"),
            "AGAL6S": Monomer(name="Gal6S", config=Monomer.Config.ALPHA,
                              smiles="O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"),
            "BGAL6S": Monomer(name="Gal6S", config=Monomer.Config.BETA,
                              smiles="O=S(=O)([O-])OC[C@@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"),

            "GALNAC": Monomer(name="GalNAc", config=Monomer.Config.UNDEF,
                              smiles="CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O"),
            "AGALNAC": Monomer(name="GalNAc", config=Monomer.Config.ALPHA,
                               smiles="CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1O"),
            "BGALNAC": Monomer(name="GalNAc", config=Monomer.Config.BETA,
                               smiles="CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1O"),

            "GALNAC4S": Monomer(name="GalNAc4S", config=Monomer.Config.UNDEF,
                                smiles="CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"),
            "AGALNAC4S": Monomer(name="GalNAc4S", config=Monomer.Config.ALPHA,
                                 smiles="CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"),
            "BGALNAC4S": Monomer(name="GalNAc4S", config=Monomer.Config.BETA,
                                 smiles="CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"),

            "GALNAC6S": Monomer(name="GalNAc6S", config=Monomer.Config.UNDEF,
                                smiles="CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"),
            "AGALNAC6S": Monomer(name="GalNAc6S", config=Monomer.Config.ALPHA,
                                 smiles="CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"),
            "BGALNAC6S": Monomer(name="GalNAc6S", config=Monomer.Config.BETA,
                                 smiles="CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"),

            "GLC": Monomer(name="Glc", config=Monomer.Config.UNDEF,
                           smiles="OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
            "AGLC": Monomer(name="Glc", config=Monomer.Config.ALPHA,
                            smiles="OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
            "BGLC": Monomer(name="Glc", config=Monomer.Config.BETA,
                            smiles="OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),

            "GLCA": Monomer(name="GlcA", config=Monomer.Config.UNDEF,
                            smiles="O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
            "AGLCA": Monomer(name="GlcA", config=Monomer.Config.ALPHA,
                             smiles="O=C(O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
            "BGLCA": Monomer(name="GlcA", config=Monomer.Config.BETA,
                             smiles="O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),

            "GLCNAC": Monomer(name="GlcNAc", config=Monomer.Config.UNDEF,
                              smiles="CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"),
            "AGLCNAC": Monomer(name="GlcNAc", config=Monomer.Config.UNDEF,
                               smiles="CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"),
            "BGLCNAC": Monomer(name="GlcNAc", config=Monomer.Config.UNDEF,
                               smiles="CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"),

            "GLCNAC6S": Monomer(name="GlcNAc6S", config=Monomer.Config.UNDEF,
                                smiles="CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"),
            "AGLCNAC6S": Monomer(name="GlcNAc6S", config=Monomer.Config.UNDEF,
                                 smiles="CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"),
            "BGLCNAC6S": Monomer(name="GlcNAc6S", config=Monomer.Config.UNDEF,
                                 smiles="CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"),

            "KDN": Monomer(name="Kdn", config=Monomer.Config.UNDEF,
                           smiles="O=C(O)C1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"),
            "AKDN": Monomer(name="Kdn", config=Monomer.Config.UNDEF,
                            smiles="O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"),
            "BKDN": Monomer(name="Kdn", config=Monomer.Config.UNDEF,
                            smiles="O=C(O)[C@]1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"),

            "MAN": Monomer(name="Man", config=Monomer.Config.UNDEF,
                           smiles="OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"),
            "AMAN": Monomer(name="Man", config=Monomer.Config.ALPHA,
                            smiles="OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"),
            "BMAN": Monomer(name="Man", config=Monomer.Config.BETA,
                            smiles="OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"),

            "NEU5AC": Monomer(name="Neu5Ac", config=Monomer.Config.UNDEF,
                              smiles="CC(=O)N[C@@H]1[C@@H](O)CC(O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"),
            "ANEU5AC": Monomer(name="Neu5Ac", config=Monomer.Config.ALPHA,
                               smiles="CC(=O)N[C@@H]1[C@@H](O)C[C@@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"),
            "BNEU5AC": Monomer(name="Neu5Ac", config=Monomer.Config.BETA,
                               smiles="CC(=O)N[C@@H]1[C@@H](O)C[C@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"),

            "NEU5GC": Monomer(name="Neu5Ac", config=Monomer.Config.UNDEF,
                              smiles="O=C(CO)N[C@@H]1[C@@H](O)CC(O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"),
            "ANEU5GC": Monomer(name="Neu5Ac", config=Monomer.Config.ALPHA,
                               smiles="O=C(CO)N[C@@H]1[C@@H](O)C[C@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"),
            "BNEU5GC": Monomer(name="Neu5Ac", config=Monomer.Config.BETA,
                               smiles="O=C(CO)N[C@@H]1[C@@H](O)C[C@@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"),

            "TAL": Monomer(name="Tal", config=Monomer.Config.UNDEF,
                           smiles="OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O"),
            "ATAL": Monomer(name="Tal", config=Monomer.Config.ALPHA,
                            smiles="OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O"),
            "BTAL": Monomer(name="Tal", config=Monomer.Config.BETA,
                            smiles="OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O"),
        }[mono.upper()]
