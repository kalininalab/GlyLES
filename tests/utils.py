import sys

from glyles.converter import convert

smiles_samples = [
    ("Gal(a1-4)Gal",
     "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),
    ("Man(a1-3)[Man(a1-6)]Man",
     "OC[C@H]3O[C@H](OC[C@H]2OC(O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H]2O)[C@@H](O)"
     "[C@@H](O)[C@@H]3O"),
    ("Man(a1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man",
     "OC[C@H]5O[C@H](O[C@H]4[C@H](O)[C@@H](CO)O[C@H](OC[C@H]3OC(O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)"
     "[C@H](O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H]3O)[C@H]4O)[C@@H](O)[C@@H](O)[C@@H]5O"),
    ("Gal(b1-4)Glc",
     "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O"),
    ("Gal(b1-3)Glc",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@H]2O"),
    ("Glc(b1-3)Glc",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O"),
    ("Gal(b1-4)Gal",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O"),
    ("Man(a1-4)Man",
     "OC[C@H]2O[C@H](O[C@H]1[C@H](O)[C@H](O)C(O)O[C@@H]1CO)[C@@H](O)[C@@H](O)[C@@H]2O"),
    ("Man(a1-3)Man",
     "OC[C@H]2O[C@H](O[C@@H]1[C@H](O)C(O)O[C@H](CO)[C@H]1O)[C@@H](O)[C@@H](O)[C@@H]2O"),
    ("Gal(a1-4)Gal(b1-4)Glc",
     "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
     "[C@@H](O)[C@H]3O"),
    ("Gal(a1-3)Gal",
     "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O)[C@H](O)[C@@H](O)[C@H]2O"),
    ("Glc(a1-3)Glc",
     "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O"),
    ("Man(a1-2)Man",
     "OC[C@H]2OC(O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)[C@@H]2O"),
]

derivatives = {
    "Gal3S": "O=S(=O)(O)O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O",
    "Gal3S a": "O=S(=O)(O)O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1O",
    "Gal3S b": "O=S(=O)(O)O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O",

    "Gal3S4S": "O=S(=O)(O)O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1OS(=O)(=O)O",
    "Gal3S4S a": "O=S(=O)(O)O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)O",
    "Gal3S4S b": "O=S(=O)(O)O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)O",

    "Gal3S6S": "O=S(=O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](OS(=O)(=O)O)[C@H]1O",
    "Gal3S6S a": "O=S(=O)(O)OC[C@H]1O[C@H](O)[C@H](O)[C@@H](OS(=O)(=O)O)[C@H]1O",
    "Gal3S6S b": "O=S(=O)(O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](OS(=O)(=O)O)[C@H]1O",

    "Gal4S": "O=S(=O)(O)O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO",
    "Gal4S a": "O=S(=O)(O)O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO",
    "Gal4S b": "O=S(=O)(O)O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO",

    "Gal4S6S": "O=S(=O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)O",
    "Gal4S6S a": "O=S(=O)(O)OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)O",
    "Gal4S6S b": "O=S(=O)(O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)O",

    "Gal6S": "O=S(=O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
    "Gal6S a": "O=S(=O)(O)OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "Gal6S b": "O=S(=O)(O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",

    "GalNAc": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O",
    "GalNAc a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1O",
    "GalNAc b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1O",

    "GalNAc3S": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)O",
    "GalNAc3S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)O",
    "GalNAc3S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)O",

    "GalNAc4S": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O",
    "GalNAc4S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O",
    "GalNAc4S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O",

    "GalNAc4S6S": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@H](OS(=O)(=O)O)[C@@H]1O",
    "GalNAc4S6S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)O)[C@H](OS(=O)(=O)O)[C@@H]1O",
    "GalNAc4S6S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)O)[C@H](OS(=O)(=O)O)[C@@H]1O",

    "GalNAc6S": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@H](O)[C@@H]1O",
    "GalNAc6S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)O)[C@H](O)[C@@H]1O",
    "GalNAc6S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)O)[C@H](O)[C@@H]1O",

    # Glucose
    "Glc4S": "O=S(=O)(O)O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO",
    "Glc4S a": "O=S(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO",
    "Glc4S b": "O=S(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO",

    "Glc6S": "O=S(=O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "Glc6S a": "O=S(=O)(O)OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "Glc6S b": "O=S(=O)(O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",

    "GlcA": "O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "GlcA a": "O=C(O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "GlcA b": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",

    "GlcN": "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
    "GlcN a": "N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",
    "GlcN b": "N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",

    "GlcNAc": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
    "GlcNAc a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",
    "GlcNAc b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",

    "GlcNAc3S": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)O",
    "GlcNAc3S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)O",
    "GlcNAc3S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)O",

    "GlcNAc6S": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",
    "GlcNAc6S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",
    "GlcNAc6S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",
}

smiles_samples_simple = [
    ("Gal(a1-4)Gal",
     "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O",
     "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O",  # alpha
     "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),  # beta

    ("Man(a1-3)[Man(a1-6)]Man",
     "OC[C@H]3O[C@H](OC[C@H]2OC(O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H]2O)[C@@H](O)"
     "[C@@H](O)[C@@H]3O",
     "OC[C@H]1O[C@H](OC[C@H]2O[C@H](O)[C@@H](O)[C@@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H]2O)[C@@H](O)"
     "[C@@H](O)[C@@H]1O",  # alpha
     "OC[C@H]1O[C@H](OC[C@H]2O[C@@H](O)[C@@H](O)[C@@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H]2O)[C@@H]"
     "(O)[C@@H](O)[C@@H]1O"),  # beta

    ("Man(a1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man",
     "OC[C@H]5O[C@H](O[C@H]4[C@H](O)[C@@H](CO)O[C@H](OC[C@H]3OC(O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H]"
     "(O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H]3O)[C@H]4O)[C@@H](O)[C@@H](O)[C@@H]5O",
     "OC[C@H]5O[C@H](O[C@H]4[C@H](O)[C@@H](CO)O[C@H](OC[C@H]3O[C@H](O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)"
     "[C@H](O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H]3O)[C@H]4O)[C@@H](O)[C@@H](O)[C@@H]5O",  # a
     "OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]2O[C@@H]2[C@H](O)[C@H](O)O[C@H](CO[C@H]3O[C@H](CO)"
     "[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),  # b

    ("Gal(b1-4)Glc",
     "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",
     "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",  # alpha
     "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O"),  # beta

    ("Gal(b1-3)Glc",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@H]2O",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@H]2O",  # alpha
     "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),  # beta

    ("Glc(b1-3)Glc",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",  # alpha
     "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),  # beta

    ("Gal(b1-4)Gal",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",
     "OC[C@H]2O[C@@H](O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",  # alpha
     "OC[C@H]1O[C@@H](O[C@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),  # beta

    ("Man(a1-4)Man",
     "OC[C@H]2O[C@H](O[C@H]1[C@H](O)[C@H](O)C(O)O[C@@H]1CO)[C@@H](O)[C@@H](O)[C@@H]2O",
     "OC[C@H]2O[C@H](O[C@H]1[C@H](O)[C@H](O)[C@@H](O)O[C@@H]1CO)[C@@H](O)[C@@H](O)[C@@H]2O",  # alpha
     "OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),  # beta

    ("Man(a1-3)Man",
     "OC[C@H]2O[C@H](O[C@@H]1[C@H](O)C(O)O[C@H](CO)[C@H]1O)[C@@H](O)[C@@H](O)[C@@H]2O",
     "OC[C@H]2O[C@H](O[C@@H]1[C@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@@H](O)[C@@H](O)[C@@H]2O",  # alpha
     "OC[C@H]1O[C@H](O[C@@H]2[C@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),  # beta

    ("Gal(a1-4)Gal(b1-4)Glc",
     "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
     "[C@@H](O)[C@H]3O",
     "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
     "[C@@H](O)[C@H]3O",  # alpha
     "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
     "[C@@H](O)[C@H]3O"),  # beta

    ("Gal(a1-3)Gal",
     "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O)[C@H](O)[C@@H](O)[C@H]2O",
     "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H]1O",  # alpha
     "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O)[C@H](O)[C@@H](O)[C@H]2O"),  # beta

    ("Glc(a1-3)Glc",
     "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",
     "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",  # alpha
     "OC[C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),  # beta

    ("Man(a1-2)Man",
     "OC[C@H]2OC(O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)[C@@H]2O",
     "OC[C@H]2O[C@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)[C@@H]2O",  # alpha
     "OC[C@H]1O[C@@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H](O)[C@@H]1O"),  # beta
]


def setup_test():
    def generator(x_list):
        for x in x_list:
            yield x

    def write_file(x_list, name="./test.txt"):
        f = open(name, "w")
        for x in x_list:
            f.write(x + "\n")
        f.close()
        return name

    return {
        "glycan": smiles_samples[0][0],
        "glycan_list": [smiles_samples[1][0], smiles_samples[2][0]],
        "glycan_file": write_file([smiles_samples[i][0] for i in range(3, 8)]),
        "glycan_generator": generator([smiles_samples[i][0] for i in range(8, 13)])
    }


def reduce_notation(iupac, mode):
    if mode != "full":
        iupac = iupac.replace("(", "").replace(")", "")
    if mode == "simple":
        iupac = iupac.replace("1-", "")
    return iupac
