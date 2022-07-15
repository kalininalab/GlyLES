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

    @pytest.mark.parametrize(
        "line",
        # open("data/glycowork_mono.txt", "r").readlines() +
        # open("data/glycowork_poly.txt", "r").readlines() +
        # open("data/general.tsv", "r").readlines() +
        # open("data/pubchem_mono.tsv", "r").readlines() +
        # open("data/pubchem_poly.tsv", "r").readlines() +
        open("data/glycam.tsv", "r").readlines()
    )
    def test_conversion_rate(self, line):
        line = line.strip()
        if '-ulosaric' in line \
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

    @pytest.mark.todo
    def test_en(self):
        iupac, smiles, _ = "Neu5Ac2en\tCC(=O)N[C@@H]1[C@H](C=C(O[C@H]1[C@@H]([C@@H](CO)O)O)C(=O)O)O\t65309".split("\t")
        output = convert(iupac)

        assert output[0][0] == iupac
        assert output[0][1] != ""
        compare_smiles(output[0][1], smiles)
