import pytest

from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdchem import ChiralType

from glyles.glycans.utils import find_isomorphism_nx as iso

CHIRALITY = [ChiralType.CHI_TETRAHEDRAL_CW, ChiralType.CHI_TETRAHEDRAL_CCW]


def check_map(smiles1, smiles2, mapping):
    mol1 = MolFromSmiles(smiles1)
    mol2 = MolFromSmiles(smiles2)

    for k, v in mapping.items():
        atom1 = mol1.GetAtomWithIdx(k)
        atom2 = mol2.GetAtomWithIdx(v)

        assert atom1.GetAtomicNum() == atom2.GetAtomicNum()


def compare_smiles(computed, solution):
    c = Chem.MolFromSmiles(computed)
    Chem.Kekulize(c)
    c_rdkit = Chem.MolToSmiles(c, kekuleSmiles=True)

    s = Chem.MolFromSmiles(solution)
    Chem.Kekulize(s)
    s_rdkit = Chem.MolToSmiles(s, kekuleSmiles=True)

    assert c_rdkit == s_rdkit


@pytest.mark.parametrize("data", [
    ("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", 12, "Glc"),  # Glc, Glc
    ("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "OC[C@H]1OC(O)C[C@@H](O)[C@@H]1O", 11, "Glc"),  # Glc, 2dGlc
    ("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "OC[C@H]1OC(O)[C@H](O)C[C@@H]1O", 11, "Glc"),  # Glc, 3dGlc
    ("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "OC[C@@H]1C[C@H](O)[C@@H](O)C(O)O1", 11, "Glc"),  # Glc, 4dGlc
    ("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "CC(=O)N[C@@H]1[C@H]([C@@H]([C@H](OC1O)CO)O)O", 11, "Glc"),  # Glc, GlcNAc
])
def test_isomorphism(data):
    smiles1, smiles2, size, name = data
    mapping = iso(smiles1, smiles2, name)

    assert len(mapping) == size
    check_map(smiles1, smiles2, mapping)
