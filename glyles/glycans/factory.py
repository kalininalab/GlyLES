from glyles.glycans.monomer import Monomer


class MonomerFactory:
    __monomers = {
        "FUC": Monomer(name="Fuc", config=Monomer.Config.UNDEF, isomer=Monomer.Enantiomer.L,
                       smiles="C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O"),
        "AFUC": Monomer(name="Fuc", config=Monomer.Config.ALPHA, isomer=Monomer.Enantiomer.L,
                        smiles="C[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O", ),
        "BFUC": Monomer(name="Fuc", config=Monomer.Config.BETA, isomer=Monomer.Enantiomer.L,
                        smiles="C[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"),

        "GAL": Monomer(name="Gal", config=Monomer.Config.UNDEF, isomer=Monomer.Enantiomer.D,
                       smiles="OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"),
        "AGAL": Monomer(name="Gal", config=Monomer.Config.ALPHA, isomer=Monomer.Enantiomer.D,
                        smiles="OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"),
        "BGAL": Monomer(name="Gal", config=Monomer.Config.BETA, isomer=Monomer.Enantiomer.D,
                        smiles="OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"),

        "GAL3S": Monomer(name="Gal3S", config=Monomer.Config.UNDEF, isomer=Monomer.Enantiomer.D,
                         smiles="O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O"),
        "AGAL3S": Monomer(name="Gal3S", config=Monomer.Config.ALPHA, isomer=Monomer.Enantiomer.D,
                          smiles="O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1O"),
        "BGAL3S": Monomer(name="Gal3S", config=Monomer.Config.BETA, isomer=Monomer.Enantiomer.D,
                          smiles="O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O"),

        "GAL4S": Monomer(name="Gal4S", config=Monomer.Config.UNDEF, isomer=Monomer.Enantiomer.D,
                         smiles="O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO"),
        "AGAL4S": Monomer(name="Gal4S", config=Monomer.Config.ALPHA, isomer=Monomer.Enantiomer.D,
                          smiles="O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO"),
        "BGAL4S": Monomer(name="Gal4S", config=Monomer.Config.BETA, isomer=Monomer.Enantiomer.D,
                          smiles="O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO"),

        "GAL6S": Monomer(name="Gal6S", config=Monomer.Config.UNDEF, isomer=Monomer.Enantiomer.D,
                         smiles="O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"),
        "AGAL6S": Monomer(name="Gal6S", config=Monomer.Config.ALPHA, isomer=Monomer.Enantiomer.D,
                          smiles="O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"),
        "BGAL6S": Monomer(name="Gal6S", config=Monomer.Config.BETA, isomer=Monomer.Enantiomer.D,
                          smiles="O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"),

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
    }

    def __contains__(self, item):
        return item.upper() in MonomerFactory.__monomers.keys()

    @staticmethod
    def __getitem__(self, item):
        return MonomerFactory.__monomers[item.upper()]
