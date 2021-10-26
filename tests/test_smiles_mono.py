import pytest

from glyles.glycans.glycans import Glycan
from glyles.smiles.smiles import SMILES

from rdkit import Chem


class TestSMILESMono:

    def compare_smiles(self, computed, solution, equal=True):
        c = Chem.MolFromSmiles(computed)
        Chem.Kekulize(c)
        c_rdkit = Chem.MolToSmiles(c, kekuleSmiles=True)

        s = Chem.MolFromSmiles(solution)
        Chem.Kekulize(s)
        s_rdkit = Chem.MolToSmiles(s, kekuleSmiles=True)

        if equal:
            assert c_rdkit == s_rdkit
        else:
            assert c_rdkit != s_rdkit

    def test_sanity(self):
        self.compare_smiles(Glycan.from_string("Man").structure(), Glycan.from_string("Glc"), equal=False)

    @pytest.mark.parameterize("glycan", ["Gly", "Fru", "Man", "Gal", "Tal"])
    @pytest.mark.parameterize("source", [1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15])
    def test_smiles_glc(self, glycan, source):
        monomer = Glycan.from_string(glycan)
        solution = monomer.value["smiles"]
        computed = SMILES().write(monomer.structure(), source)

        self.compare_smiles(computed, solution)
