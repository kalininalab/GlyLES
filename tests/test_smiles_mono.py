import pytest
from rdkit import Chem

from glyles.glycans.glycans import Glycan
from glyles.grammar.parse import parse
from glyles.smiles.smiles import SMILES, Merger


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
         "OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]2O[C@@H]2[C@H](O)[C@H](O)O[C@H](CO[C@H]3O[C@H](CO)"
         "[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O",  # a
         "OC[C@H]5O[C@H](O[C@H]4[C@H](O)[C@@H](CO)O[C@H](OC[C@H]3O[C@H](O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)"
         "[C@H](O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H]3O)[C@H]4O)[C@@H](O)[C@@H](O)[C@@H]5O"),  # b

        ("Gal(b1-4)Glc",
         "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",
         "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",  # alpha
         "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O"),  # beta

        ("Gal(b1-3)Glc",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@H]2O",
         "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O",  # alpha
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@H]2O"),  # beta

        ("Glc(b1-3)Glc",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",
         "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",  # alpha
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O"),  # beta

        ("Gal(b1-4)Gal",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O",
         "OC[C@H]1O[C@@H](O[C@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O",  # alpha
         "OC[C@H]2O[C@@H](O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O"),  # beta

        ("Man(a1-4)Man",
         "OC[C@H]2O[C@H](O[C@H]1[C@H](O)[C@H](O)C(O)O[C@@H]1CO)[C@@H](O)[C@@H](O)[C@@H]2O",
         "OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha
         "OC[C@H]2O[C@H](O[C@H]1[C@H](O)[C@H](O)[C@@H](O)O[C@@H]1CO)[C@@H](O)[C@@H](O)[C@@H]2O"),  # beta

        ("Man(a1-3)Man",
         "OC[C@H]2O[C@H](O[C@@H]1[C@H](O)C(O)O[C@H](CO)[C@H]1O)[C@@H](O)[C@@H](O)[C@@H]2O",
         "OC[C@H]1O[C@H](O[C@@H]2[C@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha
         "OC[C@H]2O[C@H](O[C@@H]1[C@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@@H](O)[C@@H](O)[C@@H]2O"),  # beta

        ("Gal(a1-4)Gal(b1-4)Glc",
         "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
         "[C@@H](O)[C@H]3O",
         "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]1O[C@H](CO)[C@H](O[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)"
         "[C@H](O)[C@H]1O",  # alpha
         "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
         "[C@@H](O)[C@H]3O"),  # beta

        ("Gal(a1-3)Gal",
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O)[C@H](O)[C@@H](O)[C@H]2O",
         "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H]1O",  # alpha
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O)[C@H](O)[C@@H](O)[C@H]2O"),  # beta

        ("Glc(a1-3)Glc",
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O",
         "OC[C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",  # alpha
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O"),  # beta

        ("Man(a1-2)Man",
         "OC[C@H]2OC(O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)[C@@H]2O",
         "OC[C@H]1O[C@@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H](O)[C@@H]1O",  # alpha
         "OC[C@H]2O[C@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)[C@@H]2O"),  # beta
    ]

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
        self.compare_smiles(Glycan.from_string("Man").value["smiles"], Glycan.from_string("Glc").value["smiles"],
                            equal=False)

    @pytest.mark.parametrize("glycan", ["Glc", "Fru", "Man", "Gal", "Tal"])
    @pytest.mark.parametrize("source", [10, 11, 12, 13, 14, 15])
    def test_smiles_mono(self, glycan, source):
        monomer = Glycan.from_string(glycan)
        solution = monomer.value["smiles"]
        computed = SMILES().write(monomer.structure(), source)

        self.compare_smiles(computed, solution)

    # @pytest.mark.parametrize("root_orientation", ["n", "a", "b"])
    @pytest.mark.parametrize("iupac, plain, alpha, beta", smiles_samples_simple)
    def test_smiles_poly(self, iupac, plain, alpha, beta, root_orientation="n"):
        computed = Merger().merge(parse(iupac, mode="nx"), root_orientation=root_orientation)

        if root_orientation == "a":
            smiles = alpha
        elif root_orientation == "b":
            smiles = beta
        else:
            smiles = plain

        self.compare_smiles(computed, smiles)
