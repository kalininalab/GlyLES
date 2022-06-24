import os
from itertools import combinations

import pytest
from rdkit import Chem

from glyles.converter import convert
from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.utils import Config
from glyles.grammar.parse import Glycan
from tests.utils import catch_output, smiles_samples_simple


def compare_smiles(computed, solution, equal=True):
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


class TestSMILES:
    def test_sanity(self):
        compare_smiles("OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
                       equal=False)

    def test_1_1_bond(self):
        smiles = Glycan("Gal(a1-1)Gal", factory=MonomerFactory()).get_smiles()
        assert smiles is not None
        assert smiles != ""

    @pytest.mark.parametrize("glycan", [
        ("Glc", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
        ("Man", "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"),
        ("Gal", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O"),
        ("Tal", "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O"),
    ])
    @pytest.mark.parametrize("source", [10, 11, 12, 13, 14, 15])
    def test_smiles_mono(self, glycan, source):
        computed = Glycan(glycan[0], factory=MonomerFactory(), start=10).get_smiles()
        solution = glycan[1]

        compare_smiles(computed, solution)

    @pytest.mark.parametrize("root_orientation", ["n", "a", "b"])
    @pytest.mark.parametrize("iupac, plain, alpha, beta", smiles_samples_simple)
    def test_smiles_poly(self, iupac, plain, alpha, beta, root_orientation):
        computed = Glycan(iupac, factory=MonomerFactory(), root_orientation=root_orientation).get_smiles()

        if root_orientation == "a":
            smiles = alpha
        elif root_orientation == "b":
            smiles = beta
        else:
            smiles = plain

        compare_smiles(computed, smiles)

    @pytest.mark.slow
    @pytest.mark.parametrize("names", list(combinations([x + "p" for x in MonomerFactory().pyranoses() if x != "Api"] +
                                                        [x + "f" for x in MonomerFactory().furanoses()], 2)))
    @pytest.mark.parametrize("config1", [Config.ALPHA, Config.BETA, Config.UNDEF])
    @pytest.mark.parametrize("config2", [Config.ALPHA, Config.BETA, Config.UNDEF])
    def test_check_different_smiles(self, names, config1, config2):
        name1, name2 = names
        factory = MonomerFactory()
        if config1 == Config.ALPHA:
            name1 += " a"
        elif config1 == Config.BETA:
            name1 += " b"

        if config2 == Config.ALPHA:
            name2 += " a"
        elif config2 == Config.BETA:
            name2 += " b"

        compare_smiles(Glycan(name1, factory).get_smiles(), Glycan(name2, factory).get_smiles(), equal=False)

    @pytest.mark.parametrize("structure", ["X?", "XX?", "?(a1-2)[X]?", "?(a1-2)[Y]?(a1-4)?", "Y?(a1-2)[X]?",
                                           "Z?(a1-2)[Y]X?"])
    def test_length(self, structure):
        arm = "?(a1-4)?(a1-4)?(a1-4)?(a1-4)?(a1-4)?(a1-4)?(a1-4)?(a1-4)?(a1-4)"
        structure = structure.replace("X", arm).replace("Y", arm[:-7]).replace("Z", arm[:-14])

        smiles = Glycan(structure.replace("?", "Gal"), MonomerFactory()).get_smiles()
        mol = Chem.MolFromSmiles(smiles)
        rings = mol.GetRingInfo().AtomRings()

        assert len(rings) == structure.count("?")
        for ring in rings:
            assert len(ring) == 6
        assert smiles.count('%') > 0
        assert smiles.count('%') % 2 == 0

        compare_smiles(smiles, Chem.MolToSmiles(mol))

    @pytest.mark.parametrize("glycan", [x[:-1] for x in open("data/glycowork_data_1.txt", "r").readlines()])
    def test_lectinoracle(self, glycan):
        if glycan in [
            'Ara(a1-5)Ara(a1-5)Ara(a1-5)Ara(a1-5)Ara(a1-5)Ara', 'Fuc(a1-3)[Neu5Ac(a2-3)]Gal(b1-4)GlcNAc(b1-3)GalNAc',
            'Fuc(a1-2)Gal(b1-3)GlcNAc(b1-3)[Fuc(a1-3)]Gal(b1-4)GlcNAc(b1-6)[Fuc(a1-2)Gal(b1-3)GlcNAc(b1-3)]'
            'Gal(b1-4)Glc', 'Ara(a1-5)Ara(a1-5)Ara(a1-5)Ara(a1-5)Ara(a1-5)Ara(a1-5)Ara',
            'Fuc(a1-2)Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-2)Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-2)Man(a1-3)]'
            'Man(b1-4)GlcNAc(b1-4)GlcNAc',
            'Gal(b1-3)GlcNAc(b1-3)[Fuc(a1-3)]Gal(b1-4)GlcNAc(b1-6)[Fuc(a1-2)Gal(b1-3)GlcNAc(b1-3)]Gal(b1-4)Glc',
            'GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)[GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)]GalNAc',
            'Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc',
            'Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-3)]GlcN(b1-3)Gal(b1-4)Glc', 'Neu5Gc(a2-3)Gal(b1-4)[Fuc(a1-4)]GlcNAc',
            'Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)]Man(a1-6)[Neu5Ac(a2-3)Gal(b1-4)'
            'GlcNAc(b1-2)Man(a1-3)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc',
        ]:
            return
        output, out, err = catch_output(convert, glycan=glycan, returning=False)

        assert output is None

        smiles = out[-2]

        mol = Chem.MolFromSmiles(smiles)
        compare_smiles(smiles, Chem.MolToSmiles(mol))

    def test_dot(self):
        glycan = Glycan("Man(a1-2)[Glc(a1-3)Gul(b1-4)]Gal(b1-3)Tal", MonomerFactory(), tree_only=True)
        glycan.save_dot("test.dot")

        with open("test.dot", "r") as dot:
            dot_lines = [x.strip() for x in dot.readlines()]

        monos = {}
        assert len([x for x in dot_lines if len(x.strip()) > 0]) == 11
        for i in range(5):
            assert dot_lines[i + 1][0].isdigit()
            assert "[label=" in dot_lines[i + 1]
            monos[dot_lines[i + 1][-5:-2]] = int(dot_lines[i + 1][0])

        assert list(sorted(monos.keys())) == ["Gal", "Glc", "Gul", "Man", "Tal"]

        edges = {
            (monos["Gal"], monos["Tal"]): "(b1-3)",
            (monos["Gul"], monos["Gal"]): "(b1-4)",
            (monos["Man"], monos["Gal"]): "(a1-2)",
            (monos["Glc"], monos["Gul"]): "(a1-3)",
        }
        for i in range(6, 10):
            assert dot_lines[i][-9:-3] == edges[(int(dot_lines[i][0]), int(dot_lines[i][5]))]

        os.remove("test.dot")

        smiles = glycan.get_smiles()
        assert smiles != ""
