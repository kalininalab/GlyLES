import os
import pytest

from glyles.converter import convert
from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.utils import ParseError
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

    @pytest.mark.parametrize("line", open("data/glycowork_mono.txt", "r").readlines())
    def test_glycowork_mono(self, line):
        iupac = line.strip()
        try:
            output = Glycan(iupac, MonomerFactory(), tree_only=True).get_tree()
        except ParseError:
            with open("data/misc/mono_not_parsed.txt", "a") as tmp:
                tmp.write(line)
                tmp.close()
            return

        assert output is not None

    @pytest.mark.parametrize("line", open("data/glycowork_poly.txt", "r").readlines())
    def test_glycowork_poly(self, line):
        iupac = line.strip()
        try:
            output = Glycan(iupac, MonomerFactory(), tree_only=True).get_tree()
        except ParseError:
            with open("data/misc/poly_not_parsed.txt", "a") as tmp:
                tmp.write(line)
                tmp.close()
            return

        assert output is not None


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
