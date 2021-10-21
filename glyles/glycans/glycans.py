from enum import Enum


class Glycan(Enum):
    GLC = {"name": "Glc", "smiles": "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O"}
    FRU = {"name": "Fru", "smiles": "C1[C@H]([C@H]([C@@H](C(O1)(CO)O)O)O)O"}
    MAN = {"name": "Man", "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"}
    GAL = {"name": "Gal", "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"}
    XYL = {"name": "Xyl", "smiles": "C1[C@H]([C@@H]([C@H](C(O1)O)O)O)O"}


def from_string(mono):
    """

    Args:
        mono (str):

    Returns:
        Glycan according to the monosaccharide provided via mono
    """
    return Glycan[mono.upper()]
