import numpy as np
from rdkit import Chem
from rdkit.Chem import GetAdjacencyMatrix, rdDepictor
from rdkit.Chem.rdchem import BondType

from glyles.glycans.utils import Config, Enantiomer, Lactole, Tree, find_longest_c_chain
from glyles.utils import smiles2mol


class OpenFactory:
    """
    Checked with http://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html
    """

    __monomers = {
        "PEN-OL": {"name": "Pent-ol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OCC(O)C(O)C(O)CO",
                   "c1_find": lambda x: c1_finder(x, "OCC(O)C(O)C(O)CO")},
        "HEX-OL": {"name": "Hex-ol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OCC(O)C(O)C(O)C(O)CO",
                   "c1_find": lambda x: c1_finder(x, "OCC(O)C(O)C(O)C(O)CO")},
        "HEP-OL": {"name": "Hept-ol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OCC(O)C(O)C(O)C(O)C(O)CO",
                   "c1_find": lambda x: c1_finder(x, "OCC(O)C(O)C(O)C(O)C(O)CO")},
        "HEPT-OL": {"name": "Hept-ol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OCC(O)C(O)C(O)C(O)C(O)CO",
                   "c1_find": lambda x: c1_finder(x, "OCC(O)C(O)C(O)C(O)C(O)CO")},
        "OCT-OL": {"name": "Oct-ol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "OCC(O)C(O)C(O)C(O)C(O)C(O)CO",
                   "c1_find": lambda x: c1_finder(x, "OCC(O)C(O)C(O)C(O)C(O)C(O)CO")},
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
        "FRU-OL": {"name": "Fructose", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)O",
                   "c1_find": lambda x: c1_finder(x, "C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)O")},
        "API-OL": {"name": "Apiol", "config": Config.UNDEF, "isomer": Enantiomer.L, "lactole": Lactole.OPEN,
                   "smiles": "OC[C@@H](O)C(O)(CO)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@@H](O)C(O)(CO)CO")},
        "ERY-OL": {"name": "Erythriol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                   "smiles": "C(O)[C@H](O)[C@H](O)CO",
                   "c1_find": lambda x: c1_finder(x, "OC[C@H](O)[C@H](O)CO")},
        "THRE-OL": {"name": "Threitol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.OPEN,
                    "smiles": "OC[C@@H](O)[C@H](O)CO",
                    "c1_find": lambda x: c1_finder(x, "OC[C@@H](O)[C@H](O)CO")},
        "INS": {"name": "Inositol", "config": Config.UNDEF, "isomer": Enantiomer.D, "lactole": Lactole.PYRANOSE,
                "smiles": "[C@H]1(O)[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1",
                "c1_find": lambda x: c1_ino_finder(x)},
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


def c1_ino_finder(structure):
    """
    For inositol, shorten that c1-finding process as this is only used for isomorphisms,
    return the list of carbon atoms.

    Args:
        structure (rdkit.Chem.rdchem.Molecule): Molecule to determine the longest carbon-chain for starting at C1

    Returns:
        List of RDKit IDs starting from C1 all the way down the longest carbon chain
    """
    a_type = np.array([a.GetAtomicNum() for a in structure.GetAtoms()])
    return list(np.where(a_type == 6)[0].tolist())


def c1_finder(structure, base_smiles):
    """
    Determine the c1 atom of monosaccharides in open form. As the "classical" method does not apply, we compare the
    structure to the SMILES string that is generated from C1.

    Args:
        structure (rdkit.Chem.rdchem.Molecule): Molecule to determine the longest carbon-chain starting at C1
        base_smiles (str): SMILES string of the root monomer starting from the oxygen of C1

    Returns:
        List of RDKit IDs starting from C1 all the way down the longest carbon chain
    """
    # compute atom types and adjacency matrix of molecule for faster computations
    a_type = np.array([a.GetAtomicNum() for a in structure.GetAtoms()])
    adjacency = GetAdjacencyMatrix(structure, useBO=True)
    c_atoms = np.where(a_type == 6)[0].tolist()

    longest_c_chain = find_longest_c_chain(c_atoms, adjacency, a_type)

    # now the two C1 candidates can be found at the ends of the longest chain
    start, end = longest_c_chain[0], longest_c_chain[-1]

    # check if exactly one of the two C1 candidates has a carboxyl group
    start_acid = check_for_acid(start, structure, adjacency, a_type)
    end_acid = check_for_acid(end, structure, adjacency, a_type)
    if start_acid and not end_acid:
        return longest_c_chain
    if end_acid and not start_acid:
        return list(reversed(longest_c_chain))

    # check if exactly one of the two C1 candidates has an attached oxygen
    start_aldehyd, end_aldehyd = check_for_aldehyd(start, adjacency, a_type), check_for_aldehyd(end, adjacency, a_type)
    if start_aldehyd and not end_aldehyd:
        return longest_c_chain
    if end_aldehyd and not start_aldehyd:
        return list(reversed(longest_c_chain))

    start_double_oxy = check_for_co_double(longest_c_chain[1], structure, adjacency, a_type)
    end_double_oxy = check_for_co_double(longest_c_chain[-2], structure, adjacency, a_type)
    if start_double_oxy and not end_double_oxy:
        return longest_c_chain
    if end_double_oxy and not start_double_oxy:
        return list(reversed(longest_c_chain))

    base = smiles2mol(base_smiles)
    if not base.GetNumConformers():
        rdDepictor.Compute2DCoords(base)
    if not structure.GetNumConformers():
        rdDepictor.Compute2DCoords(structure)
    Chem.WedgeMolBonds(base, base.GetConformer())
    Chem.WedgeMolBonds(structure, structure.GetConformer())

    base_chain = [a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() == 6]
    base_direc = [set(bond.GetBondDir() for bond in base.GetAtomWithIdx(idx).GetBonds()) for idx in base_chain]
    struct_direc = [set(bond.GetBondDir() for bond in structure.GetAtomWithIdx(idx).GetBonds()) for idx in longest_c_chain]
    equal = [bd == sd for bd, sd in zip(base_direc, struct_direc)]
    if all(equal) or not any(equal):
        return longest_c_chain
    return list(reversed(longest_c_chain))


def check_for_acid(start, structure, adjacency, a_type):
    """
    Check if the atom with the start ID has an acid group.

    Args:
        start (int): RDKit ID of the candidate for C1
        structure (rdkit.Chem.rdchem.Molecule): rdkit molecule object
        adjacency (np.ndarray): adjacency matrix of the molecule
        a_type (np.ndarray): array of the atom types of the molecule's heavy atoms

    Returns:
        bool flag indicating if the carbon atom has an acid group
    """
    oxys = [int(x) for x in np.where(np.array(adjacency[start] != 0) & (a_type == 8))[0]]
    if len(oxys) != 2:
        return False
    if (structure.GetBondBetweenAtoms(int(start), oxys[0]).GetBondType() == BondType.SINGLE and
        structure.GetBondBetweenAtoms(int(start), oxys[1]).GetBondType() == BondType.DOUBLE) or \
            (structure.GetBondBetweenAtoms(int(start), oxys[1]).GetBondType() == BondType.SINGLE and
             structure.GetBondBetweenAtoms(int(start), oxys[0]).GetBondType() == BondType.DOUBLE):
        return True
    return False


def check_for_co_double(carbon, structure, adjacency, a_type):
    oxys = [int(x) for x in np.where(np.array(adjacency[carbon] != 0) & (a_type == 8))[0]]
    return sum(structure.GetBondBetweenAtoms(int(carbon), o).GetBondType() == BondType.DOUBLE for o in oxys) == 1


def check_for_aldehyd(start, adjacency, a_type):
    """
    Check if the atom with the start ID has an aldehyd group.

    Args:
        start (int): RDKit ID of the candidate for C1
        adjacency (np.ndarray): adjacency matrix of the molecule
        a_type (np.ndarray): array of the atom types of the molecule's heavy atoms

    Returns:
        bool flag indicating if the carbon atom has an aldehyd group
    """
    return len(np.where(np.array(adjacency[start] == 1) & (a_type == 8))[0]) != 0
