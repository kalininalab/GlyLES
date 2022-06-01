from glyles.glycans.utils import *


class PyranoseFactory:
    """
    Checked with http://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html
    """

    __monomers = {
        # Hexose
        "HEX": {"name": "Glc", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"},
        "A_HEX": {"name": "Glc", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OCC1O[C@H](O)C(O)C(O)C1O"},
        "B_HEX": {"name": "Glc", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OCC1O[C@@H](O)C(O)C(O)C1O"},

        # Pentose
        "PEN": {"name": "Glc", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C1C(C(C(C(O1)O)O)O)O"},
        "A_PEN": {"name": "Glc", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC1CO[C@H](O)C(O)C1O"},
        "B_PEN": {"name": "Glc", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC1CO[C@@H](O)C(O)C1O"},

        # Glucose
        "GLC": {"name": "Glc", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_GLC": {"name": "Glc", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_GLC": {"name": "Glc", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        # Mannose
        "MAN": {"name": "Man", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "A_MAN": {"name": "Man", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "B_MAN": {"name": "Man", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},

        # Galactose
        "GAL": {"name": "Gal", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"},
        "A_GAL": {"name": "Gal", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"},
        "B_GAL": {"name": "Gal", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"},

        # Gulose
        "GUL": {"name": "Gul", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O"},
        "A_GUL": {"name": "Gul", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O"},
        "B_GUL": {"name": "Gul", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"},

        # Altrose
        "ALT": {"name": "Alt", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                "smiles": "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"},
        "A_ALT": {"name": "Alt", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"},
        "B_ALT": {"name": "Alt", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"},

        # Allose
        "ALL": {"name": "All", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O"},
        "A_ALL": {"name": "All", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@@H]1O"},
        "B_ALL": {"name": "All", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O"},

        # Talose
        "TAL": {"name": "Tal", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O"},
        "A_TAL": {"name": "Tal", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O"},
        "B_TAL": {"name": "Tal", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O"},

        # Idose
        "IDO": {"name": "Ido", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                "smiles": "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_IDO": {"name": "Ido", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_IDO": {"name": "Ido", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        # Quinovose
        "QUI": {"name": "Qui", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_QUI": {"name": "Qui", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_QUI": {"name": "Qui", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"},

        # Rhamnose
        "RHA": {"name": "Rha", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O"},
        "A_RHA": {"name": "Rha", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"},
        "B_RHA": {"name": "Rha", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O"},

        # Fucose
        "FUC": {"name": "Fuc", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O"},
        "A_FUC": {"name": "Fuc", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"},
        "B_FUC": {"name": "Fuc", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"},

        # Olivose
        "OLI": {"name": "Oli", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H]1OC(O)C[C@@H](O)[C@@H]1O"},
        "A_OLI": {"name": "Oli", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@@H](O)C[C@@H](O)[C@@H]1O"},
        "B_OLI": {"name": "Oli", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@H](O)C[C@@H](O)[C@@H]1O"},

        # Tyvelose
        "TYV": {"name": "Tyv", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H]1OC(O)[C@@H](O)C[C@@H]1O"},
        "A_TYV": {"name": "Tyv", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@H](O)[C@@H](O)C[C@@H]1O"},
        "B_TYV": {"name": "Tyv", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@@H](O)[C@@H](O)C[C@@H]1O"},

        # Abequose
        "ABE": {"name": "Abe", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H]1OC(O)[C@H](O)C[C@H]1O"},
        "A_ABE": {"name": "Abe", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@H](O)[C@H](O)C[C@H]1O"},
        "B_ABE": {"name": "Abe", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@@H](O)[C@H](O)C[C@H]1O"},

        # Paratose
        "PAR": {"name": "Par", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H]1OC(O)[C@H](O)C[C@@H]1O"},
        "A_PAR": {"name": "Par", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@H](O)[C@H](O)C[C@@H]1O"},
        "B_PAR": {"name": "Par", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@@H](O)[C@H](O)C[C@@H]1O"},

        # Digitoxose
        "DIG": {"name": "Dig", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H]1OC(O)C[C@H](O)[C@@H]1O"},
        "A_DIG": {"name": "Dig", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@H](O)C[C@H](O)[C@@H]1O"},
        "B_DIG": {"name": "Dig", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@@H](O)C[C@H](O)[C@@H]1O"},

        # Colitose
        "COL": {"name": "Col", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@@H]1OC(O)[C@@H](O)C[C@@H]1O"},
        "A_COL": {"name": "Col", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H]1O[C@H](O)[C@@H](O)C[C@@H]1O"},
        "B_COL": {"name": "Col", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H]1O[C@@H](O)[C@@H](O)C[C@@H]1O"},

        # Arabinose
        "ARA": {"name": "Ara", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                "smiles": "O[C@H]1COC(O)[C@H](O)[C@H]1O"},
        "A_ARA": {"name": "Ara", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "O[C@H]1CO[C@@H](O)[C@H](O)[C@H]1O"},
        "B_ARA": {"name": "Ara", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "O[C@H]1CO[C@H](O)[C@H](O)[C@H]1O"},

        # Lyxose
        "LYX": {"name": "Lyx", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "O[C@@H]1COC(O)[C@@H](O)[C@H]1O"},
        "A_LYX": {"name": "Lyx", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O[C@@H]1CO[C@H](O)[C@@H](O)[C@H]1O"},
        "B_LYX": {"name": "Lyx", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O[C@@H]1CO[C@@H](O)[C@@H](O)[C@H]1O"},

        # Xylose
        "XYL": {"name": "Xyl", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "O[C@@H]1COC(O)[C@H](O)[C@H]1O"},
        "A_XYL": {"name": "Xyl", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O[C@@H]1CO[C@H](O)[C@H](O)[C@H]1O"},
        "B_XYL": {"name": "Xyl", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O[C@@H]1CO[C@@H](O)[C@H](O)[C@H]1O"},

        # Ribose
        "RIB": {"name": "Rib", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "O[C@@H]1COC(O)[C@H](O)[C@@H]1O"},
        "A_RIB": {"name": "Rib", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O[C@@H]1CO[C@H](O)[C@H](O)[C@@H]1O"},
        "B_RIB": {"name": "Rib", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O[C@@H]1CO[C@@H](O)[C@H](O)[C@@H]1O"},

        # Keto-Deoxy-Nonulonic acid
        "KDN": {"name": "Kdn", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "O=C(O)C1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},
        "A_KDN": {"name": "Kdn", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},
        "B_KDN": {"name": "Kdn", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=C(O)[C@]1(O)C[C@H](O)[C@@H](O)C([C@H](O)[C@H](O)CO)O1"},

        # Neuraminic acid
        "NEU": {"name": "Neu", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "N[C@@H]1[C@@H](O)CC(O)(C(=O)O)O[C@H]1[C@H](O)[C@H](O)CO"},
        "A_NEU": {"name": "Neu", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "N[C@@H]1[C@@H](O)C[C@](O)(C(=O)O)O[C@H]1[C@H](O)[C@H](O)CO"},
        "B_NEU": {"name": "Neu", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "N[C@@H]1[C@@H](O)C[C@@](O)(C(=O)O)O[C@H]1[C@H](O)[C@H](O)CO"},

        # Pseudaminic acid
        "PSE": {"name": "Pse", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H](O)[C@H](N)C1OC(O)(C(=O)O)C[C@H](O)[C@@H]1N"},
        "A_PSE": {"name": "Pse", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H](O)[C@H](N)C1O[C@@](O)(C(=O)O)C[C@H](O)[C@@H]1N"},
        "B_PSE": {"name": "Pse", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H](O)[C@H](N)C1O[C@](O)(C(=O)O)C[C@H](O)[C@@H]1N"},

        # Legionaminic acid
        "LEG": {"name": "Leg", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@@H](O)[C@@H](N)C1OC(O)(C(=O)O)C[C@H](O)[C@H]1N"},
        "A_LEG": {"name": "Leg", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H](O)[C@@H](N)C1O[C@@](O)(C(=O)O)C[C@H](O)[C@H]1N"},
        "B_LEG": {"name": "Leg", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H](O)[C@@H](N)C1O[C@](O)(C(=O)O)C[C@H](O)[C@H]1N"},

        # Acinetaminic acid
        "ACI": {"name": "Aci", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H](O)[C@H](N)C1OC(O)(C(=O)O)C[C@H](O)[C@H]1N"},
        "A_ACI": {"name": "Aci", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H](O)[C@H](N)C1O[C@@](O)(C(=O)O)C[C@H](O)[C@H]1N"},
        "B_ACI": {"name": "Aci", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H](O)[C@H](N)C1O[C@](O)(C(=O)O)C[C@H](O)[C@H]1N"},

        # Bacillosamine
        "BAC": {"name": "Bac", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@H]1OC(O)[C@H](N)[C@@H](O)[C@@H]1N"},
        "A_BAC": {"name": "Bac", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@H](O)[C@H](N)[C@@H](O)[C@@H]1N"},
        "B_BAC": {"name": "Bac", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@H]1O[C@@H](O)[C@H](N)[C@@H](O)[C@@H]1N"},

        # L-Glycero-D-Manno-Heptose
        "LDMANHEP": {"name": "LDManHep", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "OC[C@H](O)C1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "A_LDMANHEP": {"name": "LDManHep", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "OC[C@H](O)C1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "B_LDMANHEP": {"name": "LDManHep", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "OC[C@H](O)C1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},

        # 2-Keto-3-Deoxy-D-Mannooctanoic acid
        "KDO": {"name": "Kdo", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "O=C(O)C1(O)C[C@@H](O)[C@@H](O)C([C@H](O)CO)O1"},
        "A_KDO": {"name": "Kdo", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=C(O)[C@@]1(O)C[C@@H](O)[C@@H](O)C([C@H](O)CO)O1"},
        "B_KDO": {"name": "Kdo", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=C(O)[C@@]1(O)C[C@@H](O)[C@@H](O)C([C@H](O)CO)O1"},

        # 3-Deoxy-D-Lyxo-Heptopyran-2-ularic acid
        "DHA": {"name": "Dha", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "O=C(O)[C@H]1OC(O)(C(=O)O)C[C@@H](O)[C@H]1O"},
        "A_DHA": {"name": "Dha", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=C(O)[C@H]1O[C@@](O)(C(=O)O)C[C@@H](O)[C@H]1O"},
        "B_DHA": {"name": "Dha", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "O=C(O)[C@H]1O[C@](O)(C(=O)O)C[C@@H](O)[C@H]1O"},

        # D-Glycero-D-Manno-Heptose
        "DDMANHEP": {"name": "DDManHep", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                     "smiles": "OC[C@@H](O)C1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "A_DDMANHEP": {"name": "DDManHep", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "OC[C@@H](O)C1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},
        "B_DDMANHEP": {"name": "DDManHep", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                       "smiles": "OC[C@@H](O)C1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},

        # Muramic acid
        "MUR": {"name": "Mur", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "C[C@@H](O[C@@H]1[C@@H](N)C(O)O[C@H](CO)[C@H]1O)C(=O)O"},
        "A_MUR": {"name": "Mur", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H](O[C@@H]1[C@@H](N)[C@@H](O)O[C@H](CO)[C@H]1O)C(=O)O"},
        "B_MUR": {"name": "Mur", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "C[C@@H](O[C@@H]1[C@@H](N)[C@H](O)O[C@H](CO)[C@H]1O)C(=O)O"},

        # Apiose
        "API": {"name": "Api", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.FURANOSE,
                "smiles": "OC[C@@]1(O)COC(O)[C@@H]1O"},
        "A_API": {"name": "Api", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.FURANOSE,
                  "smiles": "OC[C@@]1(O)CO[C@H](O)[C@@H]1O"},
        "B_API": {"name": "Api", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.FURANOSE,
                  "smiles": "OC[C@@]1(O)CO[C@@H](O)[C@@H]1O"},

        # Fructose
        "FRU": {"name": "Fru", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OCC1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O"},
        "A_FRU": {"name": "Fru", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O"},
        "B_FRU": {"name": "Fru", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O"},

        # Tagatose
        "TAG": {"name": "Tag", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OCC1(O)OC[C@@H](O)[C@H](O)[C@@H]1O"},
        "A_TAG": {"name": "Tag", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@]1(O)OC[C@@H](O)[C@H](O)[C@@H]1O"},
        "B_TAG": {"name": "Tag", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@@]1(O)OC[C@@H](O)[C@H](O)[C@@H]1O"},

        # Sorbose
        "SOR": {"name": "Sor", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                "smiles": "OCC1(O)OC[C@H](O)[C@@H](O)[C@@H]1O"},
        "A_SOR": {"name": "Sor", "config": Config.ALPHA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@]1(O)OC[C@H](O)[C@@H](O)[C@@H]1O"},
        "B_SOR": {"name": "Sor", "config": Config.BETA, "isomer": Enantiomer.L, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@@]1(O)OC[C@H](O)[C@@H](O)[C@@H]1O"},

        # Psicose
        "PSI": {"name": "Psi", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "OCC1(O)OC[C@@H](O)[C@@H](O)[C@H]1O"},
        "A_PSI": {"name": "Psi", "config": Config.ALPHA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@H]1O"},
        "B_PSI": {"name": "Psi", "config": Config.BETA, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                  "smiles": "OC[C@@]1(O)OC[C@@H](O)[C@@H](O)[C@H]1O"},
    }

    def __contains__(self, item):
        """
        Check if an item is part of this factory and can be returned.

        Args:
            item (str): Code of a monomer to be checked

        Returns:
            True if the item is included in the current version of this package
        """
        return item.upper() in PyranoseFactory.__monomers.keys()

    @staticmethod
    def keys():
        """
        Get all monomers that are included in this package by their extended name.

        Returns:
            Set of names for all monomers, their available derivatives and configurations (alpha/beta/undefined)
        """
        return PyranoseFactory.__monomers.keys()

    @staticmethod
    def monomers():
        """
        Get the names of all monomers in this package ignoring the alpha/beta conformations.

        Returns:
            List, sorted from long to short, of all monomer names in upper case
        """
        return list(set(x.split("_")[-1] for x in PyranoseFactory.__monomers.keys()))

    @staticmethod
    def __getitem__(item):
        """
        Get an instance of a monomer from this factory.

        Args:
            item (str): name of the query monomer

        Returns:
            Directory containing all necessary information to initialize a monomer implementation
        """
        return PyranoseFactory.__monomers[item.upper()]
