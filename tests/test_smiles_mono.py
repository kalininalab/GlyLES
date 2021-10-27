import pytest
from rdkit import Chem

from glyles.glycans.glycans import Glycan
from glyles.grammar.parse import parse
from glyles.smiles.smiles import SMILES, Merger


class TestSMILES:
    smiles_samples_simple = [
        ("Gal(a1-4)Gal",
         "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),
        # ("Gal(a1-2)Gal",  # too many Cs in SMILES! (13)
        #  "CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"),
        # ("Gal(b1-3)Gal",  # too few Cs in SMILES! (11)
        #  "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H]1O"),
        ("Man(a1-3)[Man(a1-6)]Man",
         "OC[C@H]1O[C@H](OC[C@H]2O[C@H](O)[C@@H](O)[C@@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H]2O)[C@@H](O)"
         "[C@@H](O)[C@@H]1O"),
        ("Man(a1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man",
         "OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]2O[C@@H]2[C@H](O)[C@H](O)O[C@H](CO[C@H]3O[C@H](CO)"
         "[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),
        ("Gal(b1-4)Glc",
         "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"),
        ("Gal(b1-3)Glc",
         "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),
        # ("Man(a1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man",  # too many Cs in SMILES! (32) and a C=C binding
        #  "COC[C@H]1O[C@H](OCC2=C(O)[C@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)"
        #  "[C@@H]3O)[C@H](O)[C@H](OC)O2)[C@@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H]1O"),
        ("Gal(b1-4)Gal",
         "OC[C@H]1O[C@@H](O)[C@@H]2OOC[C@H]3O[C@@H](O[C@@H]1[C@@H]2O)[C@H](O)[C@@H](O)[C@H]3O"),
        ("Glc(b1-3)Glc", "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),
        # ("Man(a1-3)[Man(a1-6)]Man",  # too many Cs in SMILES! (19) and a C=C binding?!
        #  "CO[C@H]1OC(CO[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)=C(O)[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)"
        #  "[C@@H]1O"),
        ("Gal(b1-4)Gal",
         "OC[C@H]1O[C@@H](O[C@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),
        ("Man(a1-4)Man",
         "OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),
        ("Man(a1-3)Man",
         "OC[C@H]1O[C@H](O[C@@H]2[C@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),
        ("Gal(a1-4)Gal(b1-4)Glc",
         "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]1O[C@H](CO)[C@H](O[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)"
         "[C@H](O)[C@H]1O"),
        ("Gal(a1-3)Gal",
         "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H]1O"),
        ("Glc(a1-3)Glc",
         "OC[C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),
        ("Man(a1-2)Man",
         "OC[C@H]1O[C@@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H](O)[C@@H]1O"),
        # ("Gal(a1-3)Gal(b1-4)Glc",  # Has a C=C binding
        #  "OC=C1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O[C@H]2O[C@H](CO)[C@H](O)[C@H](O)"
        #  "[C@H]2O)[C@H]1O"),
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

    @pytest.mark.parametrize("iupac, smiles", smiles_samples_simple)
    def test_smiles_poly(self, iupac, smiles):
        computed = Merger().merge(parse(iupac))

        self.compare_smiles(computed, smiles)
