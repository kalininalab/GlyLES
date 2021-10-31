from enum import Enum


class Monomer:

    class Config(Enum):
        UNDEF = 0
        ALPHA = 1
        BETA = 2

    def __init__(self):
        pass

    '''
    @staticmethod
    def from_string(mono):
        pass
    '''

    def to_smiles(self, root, ring_index):
        pass

    def alpha(self):
        pass

    def beta(self):
        pass

    def undefined(self):
        pass

    def get_config(self):
        pass

    def mark(self, position, atom):
        pass

    def get_dummy_atoms(self):
        pass

    def is_non_chiral(self):
        pass
