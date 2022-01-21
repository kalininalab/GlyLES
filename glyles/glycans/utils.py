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
    Configuration of the whole monomer regarding L and D forms in terms of enantiomerism.
    """
    D = 0
    L = 1


class Lactole(Enum):
    """
    Specification if a monomer is a pyranose (6-ring) or a furanose (5-ring)
    """
    FURANOSE = 5
    PYRANOSE = 6


class Mode(Enum):
    """
    Enumerate different modes how to represent the monomers in the tree.
    """
    DEFAULT_MODE = "rdkit"
    NETWORKX_MODE = "nx"
    RDKIT_MODE = "rdkit"


class UnreachableError(NotImplementedError):
    """
    Represent exceptions that should arise in case a piece of code is reached that under normal circumstances should
    not be reached
    """
    pass
