from enum import Enum


class Glycan(Enum):
    GLC = ("Glc")
    GLU = ("Glu")
    FRU = ("Fru")
    MAN = ("Man")
    GAL = ("Gal")


def from_string(mono: str):
    return Glycan[mono.upper()]
