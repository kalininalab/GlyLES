from glyles.glycans.utils import *


class DerivativesFactory:
    """
    Checked with http://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html
    """

    __monomers = {
        "GAL3S": {"name": "Gal3S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O"},
        "A_GAL3S": {"name": "Gal3S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1O"},
        "B_GAL3S": {"name": "Gal3S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O"},

        "GAL3S4S": {"name": "Gal3S4S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},
        "A_GAL3S4S": {"name": "Gal3S4S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                      "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},
        "B_GAL3S4S": {"name": "Gal3S4S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                      "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},

        "GAL3S6S": {"name": "Gal3S6S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},
        "A_GAL3S6S": {"name": "Gal3S6S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                      "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},
        "B_GAL3S6S": {"name": "Gal3S6S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                      "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},

        "GAL4S": {"name": "Gal4S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO"},
        "A_GAL4S": {"name": "Gal4S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO"},
        "B_GAL4S": {"name": "Gal4S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO"},

        "GAL4S6S": {"name": "Gal4S6S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},
        "A_GAL4S6S": {"name": "Gal4S6S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                      "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},
        "B_GAL4S6S": {"name": "Gal4S6S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                      "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},

        "GAL6S": {"name": "Gal6S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"},
        "A_GAL6S": {"name": "Gal6S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"},
        "B_GAL6S": {"name": "Gal6S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"},

        "GALNAC": {"name": "GalNAc", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                   "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O"},
        "A_GALNAC": {"name": "GalNAc", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1O"},
        "B_GALNAC": {"name": "GalNAc", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1O"},

        "GALNAC3S": {"name": "GalNAc3S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "A_GALNAC3S": {"name": "GalNAc3S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "B_GALNAC3S": {"name": "GalNAc3S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},

        "GALNAC4S": {"name": "GalNAc4S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},
        "A_GALNAC4S": {"name": "GalNAc4S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},
        "B_GALNAC4S": {"name": "GalNAc4S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},

        "GALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                       "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},
        "A_GALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                         "lactole": Lactole.PYRANOSE,
                         "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},
        "B_GALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.BETA, "isomer": Enantiomer.D,
                         "lactole": Lactole.PYRANOSE,
                         "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},

        "GALNAC6S": {"name": "GalNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},
        "A_GALNAC6S": {"name": "GalNAc6S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},
        "B_GALNAC6S": {"name": "GalNAc6S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},

        # Glucose
        "GLC4S": {"name": "Glc4S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO"},
        "A_GLC4S": {"name": "Glc4S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO"},
        "B_GLC4S": {"name": "Glc4S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO"},

        "GLC6S": {"name": "Glc6S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_GLC6S": {"name": "Glc6S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_GLC6S": {"name": "Glc6S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                    "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        "GLCA": {"name": "GlcA", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                 "smiles": "O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_GLCA": {"name": "GlcA", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                   "smiles": "O=C(O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_GLCA": {"name": "GlcA", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                   "smiles": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        "GLCN": {"name": "GlcN", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                 "smiles": "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "A_GLCN": {"name": "GlcN", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                   "smiles": "N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "B_GLCN": {"name": "GlcN", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                   "smiles": "N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},

        "GLCNAC": {"name": "GlcNAc", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                   "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "A_GLCNAC": {"name": "GlcNAc", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "B_GLCNAC": {"name": "GlcNAc", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},

        "GLCNAC3S": {"name": "GlcNAc3S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "A_GLCNAC3S": {"name": "GlcNAc3S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "B_GLCNAC3S": {"name": "GlcNAc3S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},

        "GLCNAC6S": {"name": "GlcNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},
        "A_GLCNAC6S": {"name": "GlcNAc6S", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},
        "B_GLCNAC6S": {"name": "GlcNAc6S", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},
    }

    def __contains__(self, item):
        """
        Check if an item is part of this factory and can be returned

        Args:
            item (str): Code of a monomer to be checked

        Returns:
            True if the item is included in the current version of this package
        """
        return item.upper() in DerivativesFactory.__monomers.keys()

    @staticmethod
    def keys():
        """
        Get all monomers that are included in this package by their extended name.

        Returns:
            Set of names for all monomer derivatives (alpha/beta/undefined)
        """
        return DerivativesFactory.__monomers.keys()

    @staticmethod
    def monomers():
        """
        Get the names of all monomers in this factory ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in upper case
        """
        return list(set(x.split("_")[-1] for x in DerivativesFactory.__monomers.keys()))

    @staticmethod
    def __getitem__(item):
        """
        Get an instance of a monomer from this factory.

        Args:
            item (str): name of the query monomer

        Returns:
            Directory containing all necessary information to initialize a monomer implementation
        """
        return DerivativesFactory.__monomers[item.upper()]
