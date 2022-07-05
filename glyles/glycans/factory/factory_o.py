from rdkit import Chem
from rdkit.Chem.rdchem import BondType

from glyles.glycans.utils import *
import numpy as np


class OpenFactory:
    """
    Checked with http://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html
    """

    __monomers = {
        "HEX-OL": {"name": "Hexol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OCC(O)C(O)C(O)C(O)CO",
                   "c1_find": lambda x: c1_finder(x, "OCC(O)C(O)C(O)C(O)CO")},
        "GLC-OL": {"name": "Glucitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO")},
        "MAN-OL": {"name": "Mannitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO")},
        "GAL-OL": {"name": "Galactitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO")},
        "GUL-OL": {"name": "Gulitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@H](O)[C@@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@H](O)[C@@H](O)[C@H](O)CO")},
        "ALT-OL": {"name": "Altritol", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@@H](O)[C@@H](O)[C@@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@@H](O)[C@@H](O)[C@@H](O)CO")},
        "ALL-OL": {"name": "Allitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@H](O)[C@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@H](O)[C@H](O)[C@H](O)CO")},
        "TAL-OL": {"name": "Talitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO")},
        "IDO-OL": {"name": "Iditol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO")},
        "QUI-OL": {"name": "Quinovitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C")},
        "RHA-OL": {"name": "Rhamnitol", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)C",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)C")},
        "FUC-OL": {"name": "Fucitol", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)C",
                   "c1_find": lambda x: c1_finder(x, "OC[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)C")},
        "ARA-OL": {"name": "Arabinitol", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@@H](O)[C@@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@@H](O)[C@@H](O)CO")},
        "LYX-OL": {"name": "Lyxitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@@H](O)[C@@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@@H](O)[C@@H](O)[C@H](O)CO")},
        "XYL-OL": {"name": "Xylitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@@H](O)[C@H](O)CO")},
        "RIB-OL": {"name": "Ribitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@H](O)[C@H](O)CO")},
        "KDO-OL": {"name": "Kdo alditol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC(=O)C(O)C[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC(=O)C(O)C[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO")},
        "MUR-OL": {"name": "Muramic alditol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](N)[C@@H](O[C@H](C)C(=O)O)[C@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](N)[C@@H](O[C@H](C)C(=O)O)[C@H](O)[C@H](O)CO")},
        "API-OL": {"name": "Apiol", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@@H](O)C(O)(CO)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@@H](O)C(O)(CO)CO")},
        "ERY-OL": {"name": "Erythriol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@H](O)CO")},
        "THRE-OL": {"name": "Threitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                    "smiles": "OC[C@@H](O)[C@H](O)CO",
                    "c1_find": lambda x: c1_finder(x, "OC[C@@H](O)[C@H](O)CO")},
    }

    def __contains__(self, item):
        """
        Check if an item is part of this factory and can be returned.

        Args:
            item (str): Code of a monomer to be checked

        Returns:
            True if the item is included in the current version of this package
        """
        return item.upper() in OpenFactory.__monomers.keys()

    @staticmethod
    def __getitem__(item):
        """
        Get an instance of a monomer from this factory.

        Args:
            item (str): name of the query monomer

        Returns:
            Directory containing all necessary information to initialize a monomer implementation
        """
        return OpenFactory.__monomers[item.upper()]


def c1_finder(structure, base_smiles):
    c_atoms = list(np.where(structure.x[:, 0] == 6)[0])
    # create a tree of all carbon atoms directly connected to the main ring of the monomer
    c_tree = Tree()
    stack = [(-1, c_atoms[0])]
    while len(stack) != 0:
        p_id, c_id = stack[-1]
        stack = stack[:-1]
        c_tree.add_node(c_id, p_id)

        children = np.argwhere((structure.adjacency[c_id, :] == 1) & (structure.x[:, 0] == 6))
        for c in children:
            if int(c) not in c_tree.nodes:
                stack.append((c_id, int(c)))

    # find the deepest node and rehang the tree to this node
    deepest_id, _ = c_tree.deepest_node()
    c_tree = c_tree.rehang_tree(deepest_id)
    longest_c_chain = c_tree.longest_chain()

    # now the two C1 candidates can be found at the ends of the longest chain
    start, end = longest_c_chain[0], longest_c_chain[-1]

    start_acid, end_acid = check_for_acid(start, structure), check_for_acid(end, structure)
    if start_acid and not end_acid:
        return longest_c_chain
    if end_acid and not start_acid:
        return reversed(longest_c_chain)

    start_aldehyd, end_aldehyd = check_for_aldehyd(start, structure), check_for_aldehyd(end, structure)
    if start_aldehyd and not end_aldehyd:
        return longest_c_chain
    if end_aldehyd and not start_aldehyd:
        return reversed(longest_c_chain)

    tmp = Chem.MolFromSmiles(base_smiles)
    tmp_chain = [a.GetIdx() for a in tmp.GetAtoms() if a.GetAtomicNum() == 6]
    for c, t in zip(longest_c_chain, tmp_chain):
        if structure.structure.GetAtomWithIdx(int(c)).GetChiralTag() != tmp.GetAtomWithIdx(int(t)).GetChiralTag():
            return reversed(longest_c_chain)
    return longest_c_chain


def check_for_acid(start, structure):
    oxys = [int(x) for x in np.where((structure.adjacency[start] != 0) & (structure.x[:, 0] == 8))[0]]
    if len(oxys) != 2:
        return False
    if (structure.structure.GetBondBetweenAtoms(int(start), oxys[0]).GetBondType() == BondType.SINGLE and
        structure.structure.GetBondBetweenAtoms(int(start), oxys[1]).GetBondType() == BondType.DOUBLE) or \
            (structure.structure.GetBondBetweenAtoms(int(start), oxys[1]).GetBondType() == BondType.SINGLE and
             structure.structure.GetBondBetweenAtoms(int(start), oxys[0]).GetBondType() == BondType.DOUBLE):
        return True
    return False


def check_for_aldehyd(start, structure):
    return np.where((structure.adjacency[start] == 1) & (structure.x[:, 0] == 8))[0] != 0
