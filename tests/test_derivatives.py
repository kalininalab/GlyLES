import pytest

from glyles.converter import convert, Glycan
from rdkit import Chem


valid_atomic_nums = [
    1, 6, 7, 8, 9, 15, 16, 17, 35, 53,
]


def compare_smiles(computed, solution):
    c = Chem.MolFromSmiles(computed)
    assert all([a.GetAtomicNum() in valid_atomic_nums for a in c.GetAtoms()])
    Chem.Kekulize(c)
    c_rdkit = Chem.MolToSmiles(c, kekuleSmiles=True)

    s = Chem.MolFromSmiles(solution)
    Chem.Kekulize(s)
    s_rdkit = Chem.MolToSmiles(s, kekuleSmiles=True)

    assert c_rdkit == s_rdkit


class TestDerivatives:
    @pytest.mark.slow
    @pytest.mark.parametrize(
        "line",
        open("data/anhydro.tsv", "r").readlines() +
        open("data/carbons.tsv", "r").readlines() +
        open("data/general.tsv", "r").readlines() +
        open("data/glycam.tsv", "r").readlines() +
        open("data/openforms.tsv", "r").readlines() +
        open("data/pubchem_mono.tsv", "r").readlines() +
        open("data/pubchem_poly.tsv", "r").readlines()
    )
    def test_smiles_databases(self, line):
        line = line.strip()
        if '0dHex' in line or 'en' in line or 'Ins' in line:
            return
        iupac, smiles = line.split("\t")[:2]
        compare_smiles(Glycan(iupac).get_smiles(), smiles)

    @pytest.mark.slow
    @pytest.mark.parametrize("line", open("data/glycowork.txt", "r").readlines())
    def test_iupac_databases(self, line):
        if '0dHex' in line or 'en' in line or 'Ins' in line:
            return
        iupac = line.strip()
        smiles = Glycan(iupac).get_smiles()
        assert smiles != ""

        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        assert all([a.GetAtomicNum() in valid_atomic_nums for a in mol.GetAtoms()])

    @pytest.mark.todo
    def test_en(self):
        iupac, smiles, _ = "Neu5Ac2en\tCC(=O)N[C@@H]1[C@H](C=C(O[C@H]1[C@@H]([C@@H](CO)O)O)C(=O)O)O\t65309".split("\t")
        output = convert(iupac)

        assert output[0][0] == iupac
        assert output[0][1] != ""
        compare_smiles(output[0][1], smiles)

    @pytest.mark.todo
    def test_ins(self):
        compare_smiles(
            Glycan("Ins").get_smiles(),
            "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
        )

    @pytest.mark.todo
    def test_ins1s6p(self):
        compare_smiles(
            Glycan("Ins1S6P").get_smiles(),
            "O=P(O)(O)O[C@@H]1[C@@H](OS(=O)(=O)O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O"
        )

    @pytest.mark.todo
    def test_ins2p4sa14gal(self):
        compare_smiles(
            Glycan("Ins2P4S(1-4)Gal").get_smiles(),
            "O=P(O)(O)O[C@@H]1[C@H](O[C@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H](OS(=O)(=O)O)[C@@H]1O"
        )
