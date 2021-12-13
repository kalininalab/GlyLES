from glyles.glycans.utils import Config, Enantiomer


class MonomerFactory:
    """
    Checked with http://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html
    """

    __monomers = {
        # Fucose
        "FUC": {"name": "Fuc", "config": Config.UNDEF, "isomer": Enantiomer.L,
                "smiles": "C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O"},
        "A_FUC": {"name": "Fuc", "config": Config.ALPHA, "isomer": Enantiomer.L,
                  "smiles": "C[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"},
        "B_FUC": {"name": "Fuc", "config": Config.BETA, "isomer": Enantiomer.L,
                  "smiles": "C[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"},

        # Galactose
        "GAL": {"name": "Gal", "config": Config.UNDEF, "isomer": Enantiomer.D,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"},
        "A_GAL": {"name": "Gal", "config": Config.ALPHA, "isomer": Enantiomer.D,
                  "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"},
        "B_GAL": {"name": "Gal", "config": Config.BETA, "isomer": Enantiomer.D,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"},

        "GAL3S": {"name": "Gal3S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O"},
        "A_GAL3S": {"name": "Gal3S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1O"},
        "B_GAL3S": {"name": "Gal3S", "config": Config.BETA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O"},

        "GAL3S4S": {"name": "Gal3S4S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},
        "A_GAL3S4S": {"name": "Gal3S4S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                      "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},
        "B_GAL3S4S": {"name": "Gal3S4S", "config": Config.BETA, "isomer": Enantiomer.D,
                      "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},

        "GAL3S6S": {"name": "Gal3S6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},
        "A_GAL3S6S": {"name": "Gal3S6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                      "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},
        "B_GAL3S6S": {"name": "Gal3S6S", "config": Config.BETA, "isomer": Enantiomer.D,
                      "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},

        "GAL4S": {"name": "Gal4S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO"},
        "A_GAL4S": {"name": "Gal4S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO"},
        "B_GAL4S": {"name": "Gal4S", "config": Config.BETA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO"},

        "GAL4S6S": {"name": "Gal4S6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},
        "A_GAL4S6S": {"name": "Gal4S6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                      "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},
        "B_GAL4S6S": {"name": "Gal4S6S", "config": Config.BETA, "isomer": Enantiomer.D,
                      "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},

        "GAL6S": {"name": "Gal6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"},
        "A_GAL6S": {"name": "Gal6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"},
        "B_GAL6S": {"name": "Gal6S", "config": Config.BETA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"},

        "GALNAC": {"name": "GalNAc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                   "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O"},
        "A_GALNAC": {"name": "GalNAc", "config": Config.ALPHA, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1O"},
        "B_GALNAC": {"name": "GalNAc", "config": Config.BETA, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1O"},

        "GALNAC3S": {"name": "GalNAc3S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "A_GALNAC3S": {"name": "GalNAc3S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "B_GALNAC3S": {"name": "GalNAc3S", "config": Config.BETA, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},

        "GALNAC4S": {"name": "GalNAc4S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},
        "A_GALNAC4S": {"name": "GalNAc4S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},
        "B_GALNAC4S": {"name": "GalNAc4S", "config": Config.BETA, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},

        "GALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},
        "A_GALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                         "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},
        "B_GALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.BETA, "isomer": Enantiomer.D,
                         "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},

        "GALNAC6S": {"name": "GalNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},
        "A_GALNAC6S": {"name": "GalNAc6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},
        "B_GALNAC6S": {"name": "GalNAc6S", "config": Config.BETA, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},

        # Glucose
        "GLC": {"name": "Glc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_GLC": {"name": "Glc", "config": Config.ALPHA, "isomer": Enantiomer.D,
                  "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_GLC": {"name": "Glc", "config": Config.BETA, "isomer": Enantiomer.D,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        "GLC4S": {"name": "Gal4S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO"},
        "A_GLC4S": {"name": "Gal4S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO"},
        "B_GLC4S": {"name": "Gal4S", "config": Config.BETA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO"},

        "GLC6S": {"name": "Glc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_GLC6S": {"name": "Glc6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_GLC6S": {"name": "Glc6S", "config": Config.BETA, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        "GLCA": {"name": "GlcA", "config": Config.UNDEF, "isomer": Enantiomer.D,
                 "smiles": "O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_GLCA": {"name": "GlcA", "config": Config.ALPHA, "isomer": Enantiomer.D,
                   "smiles": "O=C(O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_GLCA": {"name": "GlcA", "config": Config.BETA, "isomer": Enantiomer.D,
                   "smiles": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        "GLCN": {"name": "GlcN", "config": Config.UNDEF, "isomer": Enantiomer.D,
                 "smiles": "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "A_GLCN": {"name": "GlcN", "config": Config.ALPHA, "isomer": Enantiomer.D,
                   "smiles": "N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "B_GLCN": {"name": "GlcN", "config": Config.BETA, "isomer": Enantiomer.D,
                   "smiles": "N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},

        "GLCNAC": {"name": "GlcNAc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                   "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "A_GLCNAC": {"name": "GlcNAc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "B_GLCNAC": {"name": "GlcNAc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},

        "GLCNAC3S": {"name": "GalNAc3S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "A_GLCNAC3S": {"name": "GalNAc3S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "B_GLCNAC3S": {"name": "GalNAc3S", "config": Config.BETA, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},

        "GLCNAC6S": {"name": "GlcNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},
        "A_GLCNAC6S": {"name": "GlcNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},
        "B_GLCNAC6S": {"name": "GlcNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},

        # Mannose
        "MAN": {"name": "Man", "config": Config.UNDEF, "isomer": Enantiomer.D,
                "smiles": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "A_MAN": {"name": "Man", "config": Config.ALPHA, "isomer": Enantiomer.D,
                  "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "B_MAN": {"name": "Man", "config": Config.BETA, "isomer": Enantiomer.D,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},

        # Talose
        "TAL": {"name": "Tal", "config": Config.UNDEF, "isomer": Enantiomer.D,
                "smiles": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O"},
        "A_TAL": {"name": "Tal", "config": Config.ALPHA, "isomer": Enantiomer.D,
                  "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O"},
        "B_TAL": {"name": "Tal", "config": Config.BETA, "isomer": Enantiomer.D,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O"},
    }

    '''
    # Keto - Deoxy - Nonulonic acid
    "KDN": {"name": "Kdn", "config": Config.UNDEF, "smiles": "O=C(O)C1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},
    "A_KDN": {"name": "Kdn", "config": Config.UNDEF, "smiles": "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},
    "B_KDN": {"name": "Kdn", "config": Config.UNDEF, "smiles": "O=C(O)[C@]1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},
    
    # Neuraminic acid
    "NEU5AC": {"name": "Neu5Ac", "config": Config.UNDEF, "smiles": "CC(=O)N[C@@H]1[C@@H](O)CC(O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    "A_NEU5AC": {"name": "Neu5Ac", "config": Config.ALPHA, "smiles": "CC(=O)N[C@@H]1[C@@H](O)C[C@@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    "B_NEU5AC": {"name": "Neu5Ac", "config": Config.BETA, "smiles": "CC(=O)N[C@@H]1[C@@H](O)C[C@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},

    "NEU5GC": {"name": "Neu5Ac", "config": Config.UNDEF, "smiles": "O=C(CO)N[C@@H]1[C@@H](O)CC(O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    "A_NEU5GC": {"name": "Neu5Ac", "config": Config.ALPHA, "smiles": "O=C(CO)N[C@@H]1[C@@H](O)C[C@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    "B_NEU5GC": {"name": "Neu5Ac", "config": Config.BETA, "smiles": "O=C(CO)N[C@@H]1[C@@H](O)C[C@@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    '''

    def __contains__(self, item):
        return item.upper() in MonomerFactory.__monomers.keys()

    @staticmethod
    def keys():
        return MonomerFactory.__monomers.keys()

    @staticmethod
    def monomers():
        return list(set(x.split("_")[-1] for x in MonomerFactory.__monomers.keys()))

    @staticmethod
    def __getitem__(item):
        return MonomerFactory.__monomers[item.upper()]
