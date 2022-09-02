import pytest

from glyles.converter import convert
from glyles.glycans.factory.factory import MonomerFactory
from glyles.grammar.parse import Glycan
from tests.utils import derivatives
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
    @pytest.mark.parametrize("name", derivatives.keys())
    def test_basic(self, name):
        output = convert(name)[0][1]
        compare_smiles(output, derivatives[name])

    @pytest.mark.parametrize(
        "line",
        # open("data/glycowork_mono.txt", "r").readlines() +
        # open("data/glycowork_poly.txt", "r").readlines()  # +
        # open("data/general.tsv", "r").readlines() # +
        # open("data/pubchem_mono.tsv", "r").readlines()  # +
        open("data/pubchem_poly.tsv", "r").readlines()  # +
        # open("data/glycam.tsv", "r").readlines()
    )
    def test_databases(self, line):
        line = line.strip()
        if '-ulosaric' in line \
                or '-ulosonic' in line \
                or '-uronic' in line \
                or '-aric' in line \
                or '0dHex' in line \
                or 'Anhydro' in line \
                or 'en' in line \
                or 'Coum' in line \
                or 'Ins' in line:
            return
        if "\t" in line:
            (iupac, smiles), equal = line.split("\t")[:2], True
        else:
            iupac, smiles, equal = line, "", False

        computed = convert(iupac, returning=True)[0][1]
        if equal:
            compare_smiles(computed, smiles)
        else:
            assert computed != smiles
            assert all([a.GetAtomicNum() in valid_atomic_nums for a in Chem.MolFromSmiles(computed).GetAtoms()])

    def test_detail(self):
        iupac = "6dTal(a1-2)Rhaf"
        smiles = Glycan(iupac, MonomerFactory()).get_smiles()
        print(smiles)
        assert smiles != ""
        assert all([a.GetAtomicNum() in valid_atomic_nums for a in Chem.MolFromSmiles(smiles).GetAtoms()])

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
            Glycan("Ins", MonomerFactory()).get_smiles(),
            "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
        )

    @pytest.mark.todo
    def test_ins1s6p(self):
        compare_smiles(
            Glycan("Ins1S6P", MonomerFactory()).get_smiles(),
            "O=P(O)(O)O[C@@H]1[C@@H](OS(=O)(=O)O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O"
        )

    @pytest.mark.todo
    def test_ins2p4sa14gal(self):
        compare_smiles(
            Glycan("Ins2P4S(1-4)Gal", MonomerFactory()).get_smiles(),
            "O=P(O)(O)O[C@@H]1[C@H](O[C@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H](OS(=O)(=O)O)[C@@H]1O"
        )
