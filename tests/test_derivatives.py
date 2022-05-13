import pytest

from glyles.converter import convert
from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.utils import ParseError
from glyles.grammar.parse import Glycan
from tests.utils import derivatives, compare_smiles


class TestDerivatives:

    @pytest.mark.parametrize("name", derivatives.keys())
    def test_basic(self, name):
        output = convert(name, returning=True)[0][1]
        compare_smiles(output, derivatives[name])

    @pytest.mark.parametrize("line", open("data/tests.tsv", "r").readlines()[1:])
    def test_file_parsing(self, line):
        iupac, smiles = line.strip().split("\t")
        output = convert(iupac, returning=True)

        assert output[0][0] == iupac
        assert output[0][1] != ""

    @pytest.mark.parametrize("line", open("data/tests.tsv", "r").readlines()[1:])
    def test_file_correct(self, line):
        iupac, smiles = line.strip().split("\t")
        output = convert(iupac, returning=True)

        assert output[0][0] == iupac
        compare_smiles(output[0][1], smiles)

    @pytest.mark.parametrize("line", open("data/oracle.txt", "r").readlines())
    def test_oracle(self, line):
        iupac = line.strip()
        output = None
        try:
            output = Glycan(iupac, MonomerFactory(), tree_only=True).get_tree()
        except ParseError:
            with open("./still_not_parsed.txt", "a") as tmp:
                tmp.write(line)
                tmp.close()

        assert output is not None
