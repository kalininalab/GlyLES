import pytest

from glyles.converter import convert
from glyles.glycans.factory.factory import MonomerFactory
from glyles.grammar.parse import Glycan
from tests.utils import derivatives
from rdkit import Chem


def compare_smiles(computed, solution):
    c = Chem.MolFromSmiles(computed)
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

    @pytest.mark.parametrize("line", open("data/general.tsv", "r").readlines()[1:])
    def test_file_parsing(self, line):
        iupac, smiles = line.strip().split("\t")
        output = convert(iupac)

        assert output[0][0] == iupac
        assert output[0][1] != ""
        compare_smiles(output[0][1], smiles)

    @pytest.mark.parametrize(
        "line", open("data/glycowork_mono.txt", "r").readlines() + open("data/glycowork_poly.txt", "r").readlines(),
    )
    def test_glycowork_parse(self, line):
        assert Glycan(line.strip(), MonomerFactory(), tree_only=True).get_tree() is not None

    @pytest.mark.parametrize(
        "line",  # total (with exceptions on (for glycowork)    60098 / 60220   ( 99.80%)
        # open("data/glycowork_mono.txt", "r").readlines()  # |  2289 /  2391   ( 95.73%)
        open("data/glycowork_poly.txt", "r").readlines()  # |   49721 / 49731   ( 99.98%)
        # open("data/general.tsv", "r").readlines() +  # |        103 /   103   (100.00%)
        # open("data/pubchem_mono.tsv", "r").readlines() +  # |   147 /   147   (100.00%)
        # open("data/pubchem_poly.tsv", "r").readlines()  # |    7838 /  7838   (100.00%)
    )
    def test_conversion_rate(self, line):
        line = line.strip()
        if "(z" in line \
                or '-z' in line \
                or '-ulosaric' in line \
                or '-ulosonic' in line \
                or '-uronic' in line \
                or '-aric' in line \
                or 'en' in line \
                or 'Anhydro' in line \
                or 'Coum' in line \
                or 'Cer' in line \
                or '0dHex' in line \
                or 'Pau3Me7' in line \
                or 'Ins' in line \
                or 'Fuc1N4NBz7Et-ol' in line \
                or 'D-9dThrAltNon-onic' in line \
                or 'Pse5Am7Gra' in line:
            return
        if "\t" in line:
            iupac, smiles = line.split("\t")[:2]
            equal = True
        else:
            iupac, smiles, equal = line, "", False
        if equal:
            compare_smiles(Glycan(iupac, MonomerFactory()).get_smiles(), smiles)
        else:
            assert Glycan(iupac, MonomerFactory()).get_smiles() != smiles

    @pytest.mark.parametrize(
        "line", open("data/pubchem_mono.tsv", "r").readlines() + open("data/pubchem_poly.tsv", "r").readlines()
    )
    def test_pubchem(self, line):
        iupac, smiles, _ = line.strip().split("\t")
        output = convert(iupac)

        assert output[0][0] == iupac
        assert output[0][1] != ""
        compare_smiles(output[0][1], smiles)

    @pytest.mark.todo
    def test_en(self):
        iupac, smiles, _ = "Neu5Ac2en\tCC(=O)N[C@@H]1[C@H](C=C(O[C@H]1[C@@H]([C@@H](CO)O)O)C(=O)O)O\t65309".split("\t")
        output = convert(iupac)

        assert output[0][0] == iupac
        assert output[0][1] != ""
        compare_smiles(output[0][1], smiles)

    def test_quat(self):
        assert Glycan("Man(a1-4)[Man(a1-2)][Man(a1-3)][Man(a1-5)]Man(a1-4)Man", MonomerFactory(),
                      tree_only=True).get_tree() is not None
