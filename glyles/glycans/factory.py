from glyles.glycans.utils import Config, Enantiomer


class MonomerFactory:
    """
    Checked with http://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html
    """

    __monomers = {
        # Fucose
        "FUC": {"name": "Fuc", "config": Config.UNDEF, "isomer": Enantiomer.L,
                "smiles": "C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O"},
        "AFUC": {"name": "Fuc", "config": Config.ALPHA, "isomer": Enantiomer.L,
                 "smiles": "C[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"},
        "BFUC": {"name": "Fuc", "config": Config.BETA, "isomer": Enantiomer.L,
                 "smiles": "C[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"},

        # Galactose
        "GAL": {"name": "Gal", "config": Config.UNDEF, "isomer": Enantiomer.D,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"},
        "AGAL": {"name": "Gal", "config": Config.ALPHA, "isomer": Enantiomer.D,
                 "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"},
        "BGAL": {"name": "Gal", "config": Config.BETA, "isomer": Enantiomer.D,
                 "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"},

        "GAL3S": {"name": "Gal3S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O"},
        "AGAL3S": {"name": "Gal3S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1O"},
        "BGAL3S": {"name": "Gal3S", "config": Config.BETA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O"},

        "GAL3S4S": {"name": "Gal3S4S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},
        "AGAL3S4S": {"name": "Gal3S4S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                     "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},
        "BGAL3S4S": {"name": "Gal3S4S", "config": Config.BETA, "isomer": Enantiomer.D,
                     "smiles": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]"},

        "GAL3S6S": {"name": "Gal3S6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},
        "AGAL3S6S": {"name": "Gal3S6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                     "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},
        "BGAL3S6S": {"name": "Gal3S6S", "config": Config.BETA, "isomer": Enantiomer.D,
                     "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O"},

        "GAL4S": {"name": "Gal4S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO"},
        "AGAL4S": {"name": "Gal4S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO"},
        "BGAL4S": {"name": "Gal4S", "config": Config.BETA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO"},

        "GAL4S6S": {"name": "Gal4S6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                    "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},
        "AGAL4S6S": {"name": "Gal4S6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                     "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},
        "BGAL4S6S": {"name": "Gal4S6S", "config": Config.BETA, "isomer": Enantiomer.D,
                     "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]"},

        "GAL6S": {"name": "Gal6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"},
        "AGAL6S": {"name": "Gal6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"},
        "BGAL6S": {"name": "Gal6S", "config": Config.BETA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"},

        "GALNAC": {"name": "GalNAc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                   "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O"},
        "AGALNAC": {"name": "GalNAc", "config": Config.ALPHA, "isomer": Enantiomer.D,
                    "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1O"},
        "BGALNAC": {"name": "GalNAc", "config": Config.BETA, "isomer": Enantiomer.D,
                    "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1O"},

        "GALNAC3S": {"name": "GalNAc3S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "AGALNAC3S": {"name": "GalNAc3S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "BGALNAC3S": {"name": "GalNAc3S", "config": Config.BETA, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]"},

        "GALNAC4S": {"name": "GalNAc4S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},
        "AGALNAC4S": {"name": "GalNAc4S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},
        "BGALNAC4S": {"name": "GalNAc4S", "config": Config.BETA, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O"},

        "GALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                       "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},
        "AGALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                        "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},
        "BGALNAC4S6S": {"name": "GalNAc4S6S", "config": Config.BETA, "isomer": Enantiomer.D,
                        "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O"},

        "GALNAC6S": {"name": "GalNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},
        "AGALNAC6S": {"name": "GalNAc6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},
        "BGALNAC6S": {"name": "GalNAc6S", "config": Config.BETA, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O"},

        # Glucose
        "GLC": {"name": "Glc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "AGLC": {"name": "Glc", "config": Config.ALPHA, "isomer": Enantiomer.D,
                 "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "BGLC": {"name": "Glc", "config": Config.BETA, "isomer": Enantiomer.D,
                 "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        "GLC4S": {"name": "Gal4S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO"},
        "AGLC4S": {"name": "Gal4S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO"},
        "BGLC4S": {"name": "Gal4S", "config": Config.BETA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO"},

        "GLC6S": {"name": "Glc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                  "smiles": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "AGLC6S": {"name": "Glc6S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "BGLC6S": {"name": "Glc6S", "config": Config.BETA, "isomer": Enantiomer.D,
                   "smiles": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        "GLCA": {"name": "GlcA", "config": Config.UNDEF, "isomer": Enantiomer.D,
                 "smiles": "O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "AGLCA": {"name": "GlcA", "config": Config.ALPHA, "isomer": Enantiomer.D,
                  "smiles": "O=C(O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "BGLCA": {"name": "GlcA", "config": Config.BETA, "isomer": Enantiomer.D,
                  "smiles": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        "GLCN": {"name": "GlcN", "config": Config.UNDEF, "isomer": Enantiomer.D,
                 "smiles": "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "AGLCN": {"name": "GlcN", "config": Config.ALPHA, "isomer": Enantiomer.D,
                  "smiles": "N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "BGLCN": {"name": "GlcN", "config": Config.BETA, "isomer": Enantiomer.D,
                  "smiles": "N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},

        "GLCNAC": {"name": "GlcNAc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                   "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "AGLCNAC": {"name": "GlcNAc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                    "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},
        "BGLCNAC": {"name": "GlcNAc", "config": Config.UNDEF, "isomer": Enantiomer.D,
                    "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"},

        "GLCNAC3S": {"name": "GalNAc3S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "AGLCNAC3S": {"name": "GalNAc3S", "config": Config.ALPHA, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},
        "BGLCNAC3S": {"name": "GalNAc3S", "config": Config.BETA, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]"},

        "GLCNAC6S": {"name": "GlcNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                     "smiles": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},
        "AGLCNAC6S": {"name": "GlcNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},
        "BGLCNAC6S": {"name": "GlcNAc6S", "config": Config.UNDEF, "isomer": Enantiomer.D,
                      "smiles": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O"},

        # Mannose
        "MAN": {"name": "Man", "config": Config.UNDEF, "isomer": Enantiomer.D,
                "smiles": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "AMAN": {"name": "Man", "config": Config.ALPHA, "isomer": Enantiomer.D,
                 "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "BMAN": {"name": "Man", "config": Config.BETA, "isomer": Enantiomer.D,
                 "smiles": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},

        # Talose
        "TAL": {"name": "Tal", "config": Config.UNDEF, "isomer": Enantiomer.D,
                "smiles": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O"},
        "ATAL": {"name": "Tal", "config": Config.ALPHA, "isomer": Enantiomer.D,
                 "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O"},
        "BTAL": {"name": "Tal", "config": Config.BETA, "isomer": Enantiomer.D,
                 "smiles": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O"},
    }

    '''
    # Keto - Deoxy - Nonulonic acid
    "KDN": {"name": "Kdn", "config": Config.UNDEF, "smiles": "O=C(O)C1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},
    "AKDN": {"name": "Kdn", "config": Config.UNDEF, "smiles": "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},
    "BKDN": {"name": "Kdn", "config": Config.UNDEF, "smiles": "O=C(O)[C@]1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},
    
    # Neuraminic acid
    "NEU5AC": {"name": "Neu5Ac", "config": Config.UNDEF, "smiles": "CC(=O)N[C@@H]1[C@@H](O)CC(O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    "ANEU5AC": {"name": "Neu5Ac", "config": Config.ALPHA, "smiles": "CC(=O)N[C@@H]1[C@@H](O)C[C@@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    "BNEU5AC": {"name": "Neu5Ac", "config": Config.BETA, "smiles": "CC(=O)N[C@@H]1[C@@H](O)C[C@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},

    "NEU5GC": {"name": "Neu5Ac", "config": Config.UNDEF, "smiles": "O=C(CO)N[C@@H]1[C@@H](O)CC(O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    "ANEU5GC": {"name": "Neu5Ac", "config": Config.ALPHA, "smiles": "O=C(CO)N[C@@H]1[C@@H](O)C[C@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    "BNEU5GC": {"name": "Neu5Ac", "config": Config.BETA, "smiles": "O=C(CO)N[C@@H]1[C@@H](O)C[C@@](O)(C(=O)O)OC1[C@H](O)[C@H](O)CO"},
    '''

    def __contains__(self, item):
        return item.upper() in MonomerFactory.__monomers.keys()

    @staticmethod
    def keys():
        return MonomerFactory.__monomers.keys()

    @staticmethod
    def __getitem__(item):
        return MonomerFactory.__monomers[item.upper()]
