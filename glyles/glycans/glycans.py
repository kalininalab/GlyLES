from enum import Enum


class Glycan(Enum):
    GLC = ("Glc", "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O")
    FRU = ("Fru", "C1[C@H]([C@H]([C@@H](C(O1)(CO)O)O)O)O")
    MAN = ("Man", "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O")
    GAL = ("Gal", "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    XYL = ("Xyl", "C1[C@H]([C@@H]([C@H](C(O1)O)O)O)O")

    def __init__(self, name, smiles):
        self.name = name
        self.smiles = smiles


def from_string(mono: str):
    return Glycan[mono.upper()]
