import logging

import numpy as np
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdmolops import AddHs, RemoveHs

from glyles.glycans.utils import Enantiomer, ketoses2
from glyles.grammar.GlycanLexer import GlycanLexer


functional_groups = {
    "": "",
    "F": "F",
    "I": "I",
    "N": "N",
    "S": "S(=O)(=O)(O)",
    "P": "P(=O)(O)(O)",
    "Ac": "C(C)(=O)",
    "Cl": "Cl",
    "Bn": "c2ccccc2",
    "Br": "Br",
    "Bz": "C(=O)c2ccccc2",
    "Et": "CC",
    "Gc": "C(=O)CO",
    "Me": "C",
    "Ph": "c2ccccc2",
    "Tf": "S(=O)(=O)C(F)(F)F",
    "Tr": "C(c2ccccc2)(c3ccccc3)c4ccccc4",
    "Ts": "S(=O)(=O)c2ccc(C)cc2",
    # "Ala": "C(=O)[C@@H](C)N",
    "Asp": "N[C@@H](CC(=O)O)C(=O)O",
    "But": "CCCC",
    # "Cho": "P(=O)(O)OCC[N+](C)(C)C",
    "Cys": "N[C@@H](CS)C(=O)O",
    "Gly": "OCC(O)CO",
    "Lys": "NCCCC[C@H](N)C(=O)O",
    "Mal": "C(C(=O)O)CC(=O)O",
    "Ole": "OC(=O)CCCCCCC/C=C\CCCCCCCC",
    "Pro": "N1CCCC1C(=O)O",
    "Thr": "N[C@H](C(=O)O)[C@@H](O)C",
    "Prop": "CCC",
}


def not_implemented_message(mod):
    """
    Print a message if some modification is parsable but not implemented yet. So, the grammar might accept some
    modifications but the reactor is yet not capable to attach that functional group to the monomer

    Args
        mod (str): name of the modification that cannot be converted

    Returns
        Nothing
    """
    logging.warning(
        f"ModificationNotImplementedWarning: {mod} Modification not implemented. The returned molecule will not have "
        f"this modification")


def opposite_chirality(tag):
    if tag == ChiralType.CHI_TETRAHEDRAL_CCW:
        return ChiralType.CHI_TETRAHEDRAL_CW
    elif tag == ChiralType.CHI_TETRAHEDRAL_CW:
        return ChiralType.CHI_TETRAHEDRAL_CCW
    return tag


class SMILESReaktor:
    def __init__(self, monomer):
        self.monomer = monomer
        self.monomer.get_structure()
        self.side_chains = None
        self.ring_c = 2 if (monomer, monomer.get_lactole()) in ketoses2 else 1

    def react(self, names, types):
        full = True

        # parse for sth. like LDManHep or DDManHep -> afterwards, no need to look for LD/DD/... and Hep/Hex/Pen/Oct
        self.check_for_resizing(names, types)

        self.side_chains = ["" for _ in range(1 + np.count_nonzero(self.monomer.x[:, 0] == 6))]

        # parse remaining modifications
        for n, t in zip(names, types):
            if t != GlycanLexer.MOD or n.count("L") + n.count("D") == len(n) or n == '-':
                continue
            if n == "A":
                self.side_chains[-1] = "(=O)O"
            elif n == "N":
                self.side_chains[1 if self.monomer.get_name() in ['Fru', 'Tag', 'Sor', 'Psi'] else 2] = "N"
            elif n == "D-":
                self.to_enantiomer(Enantiomer.D)
            elif n == "L-":
                self.to_enantiomer(Enantiomer.L)
            elif n == "Ac" and self.monomer.get_name() == "Neu":
                self.side_chains[5] = "NC(C)(=O)"
            elif n == "Gc" and self.monomer.get_name() == "Neu":
                self.side_chains[5] = "NC(=O)CO"
            elif n[0].isdigit():
                if n[1:] == "d":
                    self.side_chains[int(n[0])] = "H"
                elif n[1:] == "e":  # change chirality at a single chiral carbon atom
                    idx = int(np.where(self.monomer.x[:, 1] == int(n[0]))[0])
                    tag = opposite_chirality(self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag())
                    self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(tag)
                elif n[1:4] == "-O-":
                    full &= self.set_fg(int(n[0]), "O", n[4:-1])
                elif n[1] in "ON":
                    full &= self.set_fg(int(n[0]), n[1], n[2:])
                else:
                    if n[1:] in ["Ac", "Ala", "Bz", "Bn", "But", "Cho", "Et", "Gc", "Me", "Prop", "Tf", "Ts"]:
                        elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol()
                        full &= self.set_fg(int(n[0]), elem, n[1:])
                    elif n[1:] in ["P", "S"]:
                        full &= self.set_fg(int(n[0]), "O", n[1:])
                    else:
                        full &= self.set_fg(int(n[0]), "", n[1:])
            elif n[0] in "NO" and n not in ["Ole"]:
                if n[1:] == "Me":
                    full &= self.set_fg(self.ring_c, n[0], n[1:])
                else:
                    full &= self.set_fg(self.ring_c + 1, n[0], n[1:])
            else:
                if n in ["But", "Et", "Prop"]:
                    elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(self.ring_c)).GetSymbol()
                else:
                    elem = ""
                full &= self.set_fg(self.ring_c, elem, n)

        self.assemble_chains()

        return self.monomer, full

    def set_fg(self, pos, bond_elem, name):
        if name in functional_groups:
            self.side_chains[pos] = bond_elem + functional_groups[name]
            return True
        else:
            not_implemented_message(name)
            return False

    def check_for_resizing(self, names, types):
        sac_index = types.index(GlycanLexer.SAC)
        if len(types) == sac_index + 1 or (len(types) > sac_index + 1 and types[sac_index + 1] != GlycanLexer.SAC):
            return
        if names[sac_index + 1] in ["Pen", "Hex", "Hep", "Oct"]:
            len_index = sac_index + 1
        else:
            len_index = sac_index
            sac_index = len_index + 1
        if sac_index == 0:
            orientations = ""
        else:
            orientations = names[sac_index - 1]

        if names[len_index] == "Hep":
            extension = "[C?H](O)CO"
        elif names[len_index] == "Oct":
            return
        else:
            return

        for c in orientations[:-1]:
            if c == "L":
                extension = extension.replace("[C?H]", "[C@@H]", 1)
            if c == "D":
                extension = extension.replace("[C?H]", "[C@H]", 1)

        self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(6)).SetAtomicNum(50)
        self.monomer.structure.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == 6)[0])).SetAtomicNum(32)
        smiles = self.monomer.to_smiles(int(np.where(self.monomer.x[:, 1] == 10)[0]), 0)
        smiles = smiles.replace("[SnH]", "").replace("[GeH2]", extension)
        self.monomer.smiles = smiles
        self.monomer.structure = None
        self.monomer.get_structure()

    def assemble_chains(self):
        placeholder = [
            (31, "[GaH2]"), (32, "[GeH3]"),
            (49, "[InH2]"), (50, "[SnH]"), (51, "[SbH2]"), (52, "[TeH]"),
            (81, "[TlH2]"), (82, "[PbH]"), (83, "[BiH2]"), (84, "[PoH]"),
        ]
        for i, chain in enumerate(self.side_chains):
            if chain:
                try:
                    idx = self.monomer.find_oxygen(i)
                    if chain == "H":  # account for desoxygenation
                        self.monomer.structure.GetAtomWithIdx(idx).SetAtomicNum(1)
                    else:
                        self.monomer.structure.GetAtomWithIdx(idx).SetAtomicNum(placeholder[i][0])
                except ValueError:
                    tmp = AddHs(self.monomer.structure)
                    h = None
                    for n in tmp.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == i)[0])).GetNeighbors():
                        if n.GetAtomicNum() == 1:
                            h = n.GetIdx()
                            break
                    if h is None:
                        raise ValueError(f"There is no oxygen, nitrogen, and hydrogen attached to C{i}! "
                                         f"No functional group can be attached there!")
                    tmp.GetAtomWithIdx(h).SetAtomicNum(placeholder[i][0])
                    self.monomer.structure = RemoveHs(tmp)
                    if chain[0] == "O":
                        self.side_chains[i] = self.side_chains[i][1:]
        smiles = self.monomer.to_smiles(int(np.where(self.monomer.x[:, 1] == 10)[0]), 0)
        for i, chain in enumerate(self.side_chains):
            if chain:
                smiles = smiles.replace(placeholder[i][1], "" if chain == "H" else chain)
        self.monomer.smiles = smiles.replace("()", "")
        self.monomer.structure = None
        self.monomer.get_structure()

    def to_enantiomer(self, form):
        if self.monomer.get_isomer() == form:
            return
        for idx in [atom.GetIdx() for atom in self.monomer.structure.GetAtoms()]:
            """self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(opposite_chirality(
                self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag()
            ))"""
            if self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag() == ChiralType.CHI_TETRAHEDRAL_CCW:
                self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CW)
            elif self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag() == ChiralType.CHI_TETRAHEDRAL_CW:
                self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CCW)

