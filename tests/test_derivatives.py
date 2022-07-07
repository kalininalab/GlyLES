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
        "line", open("data/glycowork_mono.txt", "r").readlines() + open("data/misc/glycowork_poly.txt", "r").readlines(),
    )
    def test_glycowork_parse(self, line):
        assert Glycan(line.strip(), MonomerFactory(), tree_only=True).get_tree() is not None

    @pytest.mark.parametrize(
        "line",  # total (with exceptions on (for glycowork)    40188 / 40198   ( 99.98%)
        # open("data/glycowork_mono.txt", "r").readlines()  # |  2391 /  2391   (100.00%)
        open("data/glycowork_poly.txt", "r").readlines()  # |   29709 / 29719   ( 99.97%)
        # open("data/general.tsv", "r").readlines() +  # |        103 /   103   (100.00%)
        # open("data/pubchem_mono.tsv", "r").readlines() +  # |   147 /   147   (100.00%)
        # open("data/pubchem_poly.tsv", "r").readlines()  # |    7838 /  7838   (100.00%)
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

    @pytest.mark.parametrize("iupac", [
        (0, "Araf(a1-3)[Gal(b1-6)]Gal(b1-5)Araf(a1-5)Araf(a1-5)Araf(a1-5)Araf(a1-5)Araf(a1-5)Araf(a1-5)Araf(a1-5)Araf"
            "(a1-6)Gal(b1-6)Gal(b1-6)Gal(b1-6)Gal(b1-6)Gal(b1-6)Gal(b1-6)Gal(b1-6)Gal(b1-6)Gal(b1-6)[Araf(a1-5)Araf"
            "(a1-5)Araf(a1-5)Araf(a1-5)Araf(a1-5)Araf(a1-3)]Gal][GlcA(a1-2)Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)]"
            "Xyl(b1-4)]Xyl(b1-4)Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)]Xyl"
            "(b1-4)]Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)]Xyl(b1-4)]Xyl"),  # invalid >]Gal][<
        (1, "DDManHep(a1-7)LDManHep(a1-3)LDManHep(a1-7)[Rha(a1-2)]LDManHep(a1-3)[Glc(b1-4)]LDManHep(a1-5)Kdo(a2-4)Kdo"
            "(b1-8)Ara4N"),  # invalid
        (2, "Gal(b1-6)Gal(b1-5)Araf(a1-6)[Gal(b1-3)Rha(a1-2)[Gal(b1-3)GalA(b1-4)][Gal(b1-6)]Araf(a1-3)]Gal(b1-3)Gal"
            "(b1-6)Gal(b1-5)Araf(a1-6)[Gal(b1-3)GlcA(b1-3)Araf(a1-3)[Gal(b1-3)Rha(a1-2)][Gal(b1-4)]Araf(a1-3)]Gal"
            "(b1-3)Gal"),  # invalid two bonds to C6 in first line
        (3, "Glc(b1-3)GalNAc(a1-3)GalNAc(b1-3)Rha(a1-3)[QuiNAc(b1-3)Rha(a1-2)][LDManHep(a1-7)]LDManHep(a1-3)"
            "[LDManHep(a1-6)Glc(b1-4)]LDManHep(a1-5)Kdo(a2-4)Kdo(b1-8)AraN"),  # invalid >-8)Ara<
        (4, "Rha(a1-3)QuiNAc(b1-7)LDManHep(a1-2)[Glc(a1-6)]Gal(a1-3)[LDManHep(a1-7)LDManHep(a1-7)]LDManHep(a1-3)"
            "[Glc(b1-4)]LDManHep(a1-5)Kdo(a2-4)Kdo(b1-8)AraN"),  # invalid >-8)Ara<
        (5, "Suc(a1-6)Glc(b1-2)GlcOPGro"),
        (6, "Suc(a1-6)Glc(b1-2)GlcPGro"),
        (7, "Gal(a1-2)GlcA(b1-4)[XylOMe(a1-3)]Fuc(a1-4)[GalA(a1-2)][GalA(b1-3)]Rha(b1-5)Apif"),
        (8, "Rha(a1-2)Ara(a1-4)[FucOMe(a1-2)]Gal(a1-2)AcefA(b1-3)Rha(b1-5)Apif"),
        (9, "XylOMe(a1-3)FucOMe(a1-4)RhaOMe(b1-5)ApiOMe-ol"),
    ])
    def test_glycowork_poly_detail(self, iupac):
        pass

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
