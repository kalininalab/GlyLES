from enum import Enum


class Config(Enum):
    """
    Configuration enum to represent if the monomer is in alpha, beta, or undefined configuration.
    """
    UNDEF = 0
    ALPHA = 1
    BETA = 2


class Enantiomer(Enum):
    """
    Configuration of the whole monomer regarding L and D forms in term of enantiomerism.
    """
    D = 0
    L = 1
