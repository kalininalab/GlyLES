import pytest

from glyles.converter import convert
from tests.test_smiles import compare_smiles


class TestDerivatives:
    __derivatives = {
        "Gal3S": "O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O",
        "Gal3S a": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1O",
        "Gal3S b": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1O",

        "Gal3S4S": "O=S(=O)([O-])O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]",
        "Gal3S4S a": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]",
        "Gal3S4S b": "O=S(=O)([O-])O[C@@H]1[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]1OS(=O)(=O)[O-]",

        "Gal3S6S": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O",
        "Gal3S6S a": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O",
        "Gal3S6S b": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@H]1O",

        "Gal4S": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO",
        "Gal4S a": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO",
        "Gal4S b": "O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO",

        "Gal4S6S": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]",
        "Gal4S6S a": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]",
        "Gal4S6S b": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1OS(=O)(=O)[O-]",

        "Gal6S": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
        "Gal6S a": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
        "Gal6S b": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",

        "GalNAc": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O",
        "GalNAc a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1O",
        "GalNAc b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1O",

        "GalNAc3S": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]",
        "GalNAc3S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]",
        "GalNAc3S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1OS(=O)(=O)[O-]",

        "GalNAc4S": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O",
        "GalNAc4S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O",
        "GalNAc4S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](OS(=O)(=O)O)[C@@H]1O",

        "GalNAc4S6S": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O",
        "GalNAc4S6S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O",
        "GalNAc4S6S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H]1O",

        "GalNAc6S": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O",
        "GalNAc6S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O",
        "GalNAc6S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@H](O)[C@@H]1O",

        # Glucose
        "Glc4S": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO",
        "Glc4S a": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1CO",
        "Glc4S b": "O=S(=O)([O-])O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]1CO",

        "Glc6S": "O=S(=O)([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Glc6S a": "O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Glc6S b": "O=S(=O)([O-])OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",

        "GlcA": "O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "GlcA a": "O=C(O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "GlcA b": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",

        "GlcN": "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
        "GlcN a": "N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",
        "GlcN b": "N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",

        "GlcNAc": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
        "GlcNAc a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",
        "GlcNAc b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",

        "GlcNAc3S": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]",
        "GlcNAc3S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]",
        "GlcNAc3S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1OS(=O)(=O)[O-]",

        "GlcNAc6S": "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O",
        "GlcNAc6S a": "CC(=O)N[C@H]1[C@@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O",
        "GlcNAc6S b": "CC(=O)N[C@H]1[C@H](O)O[C@H](COS(=O)(=O)[O-])[C@@H](O)[C@@H]1O",
    }

    @pytest.mark.parametrize("name", __derivatives.keys())
    def test_basic(self, name):
        output = convert(name, returning=True)[0][1]
        compare_smiles(output, self.__derivatives[name])

    def test_basic_detail(self, name="Gal3S"):
        output = convert(name, returning=True)[0][1]
        compare_smiles(output, self.__derivatives[name])
