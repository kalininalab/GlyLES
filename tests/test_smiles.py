from itertools import combinations

import pytest
from rdkit import Chem

from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.utils import Mode, Config
from glyles.grammar.parse import Glycan


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
    smiles_samples_simple = [
        ("Gal(a1-4)Gal",
         "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O",
         "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O",  # alpha
         "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),  # beta

        ("Man(a1-3)[Man(a1-6)]Man",
         "OC[C@H]3O[C@H](OC[C@H]2OC(O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H]2O)[C@@H](O)"
         "[C@@H](O)[C@@H]3O",
         "OC[C@H]1O[C@H](OC[C@H]2O[C@H](O)[C@@H](O)[C@@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H]2O)[C@@H](O)"
         "[C@@H](O)[C@@H]1O",  # alpha
         "OC[C@H]1O[C@H](OC[C@H]2O[C@@H](O)[C@@H](O)[C@@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H]2O)[C@@H]"
         "(O)[C@@H](O)[C@@H]1O"),  # beta

        ("Man(a1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man",
         "OC[C@H]5O[C@H](O[C@H]4[C@H](O)[C@@H](CO)O[C@H](OC[C@H]3OC(O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H]"
         "(O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H]3O)[C@H]4O)[C@@H](O)[C@@H](O)[C@@H]5O",
         "OC[C@H]5O[C@H](O[C@H]4[C@H](O)[C@@H](CO)O[C@H](OC[C@H]3O[C@H](O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)"
         "[C@H](O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H]3O)[C@H]4O)[C@@H](O)[C@@H](O)[C@@H]5O",  # a
         "OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]2O[C@@H]2[C@H](O)[C@H](O)O[C@H](CO[C@H]3O[C@H](CO)"
         "[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),  # b

        ("Gal(b1-4)Glc",
         "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",
         "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",  # alpha
         "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O"),  # beta

        ("Gal(b1-3)Glc",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@H]2O",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@H]2O",  # alpha
         "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),  # beta

        ("Glc(b1-3)Glc",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",  # alpha
         "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),  # beta

        ("Gal(b1-4)Gal",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",  # alpha
         "OC[C@H]1O[C@@H](O[C@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),  # beta

        ("Man(a1-4)Man",
         "OC[C@H]2O[C@H](O[C@H]1[C@H](O)[C@H](O)C(O)O[C@@H]1CO)[C@@H](O)[C@@H](O)[C@@H]2O",
         "OC[C@H]2O[C@H](O[C@H]1[C@H](O)[C@H](O)[C@@H](O)O[C@@H]1CO)[C@@H](O)[C@@H](O)[C@@H]2O",  # alpha
         "OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),  # beta

        ("Man(a1-3)Man",
         "OC[C@H]2O[C@H](O[C@@H]1[C@H](O)C(O)O[C@H](CO)[C@H]1O)[C@@H](O)[C@@H](O)[C@@H]2O",
         "OC[C@H]2O[C@H](O[C@@H]1[C@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@@H](O)[C@@H](O)[C@@H]2O",  # alpha
         "OC[C@H]1O[C@H](O[C@@H]2[C@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),  # beta

        ("Gal(a1-4)Gal(b1-4)Glc",
         "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
         "[C@@H](O)[C@H]3O",
         "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
         "[C@@H](O)[C@H]3O",  # alpha
         "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
         "[C@@H](O)[C@H]3O"),  # beta

        ("Gal(a1-3)Gal",
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O)[C@H](O)[C@@H](O)[C@H]2O",
         "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H]1O",  # alpha
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O)[C@H](O)[C@@H](O)[C@H]2O"),  # beta

        ("Glc(a1-3)Glc",
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",  # alpha
         "OC[C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),  # beta

        ("Man(a1-2)Man",
         "OC[C@H]2OC(O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)[C@@H]2O",
         "OC[C@H]2O[C@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)[C@@H]2O",  # alpha
         "OC[C@H]1O[C@@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H](O)[C@@H]1O"),  # beta
    ]

    def test_sanity(self):
        compare_smiles("OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
                       equal=False)

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
        computed = Glycan(iupac, factory=MonomerFactory(), mode=Mode.RDKIT_MODE,
                          root_orientation=root_orientation).get_smiles()

        if root_orientation == "a":
            smiles = alpha
        elif root_orientation == "b":
            smiles = beta
        else:
            smiles = plain

        compare_smiles(computed, smiles)

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

    @pytest.mark.parametrize("structure", ["X?", "XX?", "?(a1-2)[Y]?", "?(a1-2)[Z]?(a1-4)?", "Z?(a1-2)[Y]?",
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

        compare_smiles(smiles, Chem.MolToSmiles(mol))
