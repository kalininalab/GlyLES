import logging

import numpy as np
from rdkit.Chem import MolToSmiles
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdmolops import AddHs, RemoveHs

from glyles.glycans.utils import Enantiomer, ketoses2
from glyles.grammar.GlycanLexer import GlycanLexer


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
            if t != GlycanLexer.MOD or n.count("L") + n.count("D") == len(n):
                continue
            if len(n) == 1:
                if n == "A":
                    self.side_chains[-1] = "(=O)O"
                elif n == "N":
                    self.side_chains[1 if self.monomer.get_name() in ['Fru', 'Tag', 'Sor', 'Psi'] else 2] = "N"
                else:
                    not_implemented_message(n)
                    full = False
            elif len(n) == 2:
                if n[0].isdigit():
                    if n[1] == "d":
                        self.side_chains[int(n[0])] = "H"
                    elif n[1] == "e":  # change chirality at a single chiral carbon atom
                        idx = int(np.where(self.monomer.x[:, 1] == int(n[0]))[0])
                        tag = opposite_chirality(self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag())
                        self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(tag)
                    elif n[1] == "F":
                        self.side_chains[int(n[0])] = "F"
                    elif n[1] == "I":
                        self.side_chains[int(n[0])] = "I"
                    elif n[1] == "N":  # add a nitrogen atom to a certain position (TBT with O)
                        self.side_chains[int(n[0])] = "N"
                    elif n[1] == "P":  # add a phosphate atom to a certain position (TBT with O)
                        self.side_chains[int(n[0])] = "OP(=O)(O)(O)"
                    elif n[1] == "S":  # add a sulfur atom to a certain position (TBT with O)
                        self.side_chains[int(n[0])] = "OS(=O)(=O)(O)"
                    else:
                        not_implemented_message(n)
                        full = False
                elif n[0] in "ON":
                    if n[1] == "F":  # add a nitrogen atom to a certain position (TBT with O)
                        self.side_chains[self.ring_c + 1] = n[0] + "F"
                    if n[1] == "I":  # add a nitrogen atom to a certain position (TBT with O)
                        self.side_chains[self.ring_c + 1] = n[0] + "I"
                    if n[1] == "N":  # add a nitrogen atom to a certain position (TBT with O)
                        self.side_chains[self.ring_c + 1] = n[0] + "N"
                    elif n[1] == "S":  # add a sulfur atom to a certain position (TBT with O)
                        self.side_chains[self.ring_c + 1] = n[0] + "S(=O)(=O)(O)"
                    elif n[1] == "P":  # add a phosphate atom to a certain position (TBT with O)
                        self.side_chains[self.ring_c + 1] = n[0] + "P(=O)(O)(O)"
                    else:
                        not_implemented_message(n)
                        full = False
                elif n == "D-":
                    self.to_enantiomer(Enantiomer.D)
                elif n == "L-":
                    self.to_enantiomer(Enantiomer.L)
                elif n == "Ac" and self.monomer.get_name() == "Neu":
                    self.side_chains[5] = "NC(C)(=O)"
                elif n == "Gc" and self.monomer.get_name() == "Neu":
                    self.side_chains[5] = "NC(=O)CO"
                else:
                    not_implemented_message(n)
                    full = False
            elif len(n) == 3:
                if n[0].isdigit():
                    if n[1:] == "Me":
                        self.side_chains[int(n[0])] = "OC"
                    elif n[1:] == "Ac":
                        elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol()
                        self.side_chains[int(n[0])] = elem + "C(C)(=O)"
                    elif n[1:] == "Cl":
                        self.side_chains[int(n[0])] = "Cl"
                    elif n[1:] == "Bn":
                        self.side_chains[int(n[0])] = "Oc2ccccc2"
                    elif n[1:] == "Br":
                        self.side_chains[int(n[0])] = "Br"
                    elif n[1:] == "Bz":
                        self.side_chains[int(n[0])] = "OC(=O)c2ccccc2"
                    elif n[1:] == "Gc":
                        self.side_chains[int(n[0])] = ("N" if self.monomer.get_name() == "Neu" else "O") + "C(=O)CO"
                    elif n[1:] == "Ph":
                        self.side_chains[int(n[0])] = "c2ccccc2"
                    elif n[1:] == "Tf":
                        self.side_chains[int(n[0])] = "OS(=O)(=O)C(F)(F)F"
                    elif n[1:] == "Tr":
                        self.side_chains[int(n[0])] = "C(c2ccccc2)(c3ccccc3)c4ccccc4"
                    elif n[1:] == "Ts":
                        self.side_chains[int(n[0])] = "OS(=O)(=O)c2ccc(C)cc2"
                    else:
                        not_implemented_message(n)
                        full = False
                elif n[0] in "ON":
                    if n[1:] == "Ac":
                        self.side_chains[self.ring_c + 1] = n[0] + "C(C)(=O)"
                    elif n[1:] == "Cl":
                        self.side_chains[self.ring_c + 1] = n[0] + "Cl"
                    elif n[1:] == "Bn":
                        self.side_chains[self.ring_c + 1] = n[0] + "c2ccccc2"
                    elif n[1:] == "Br":
                        self.side_chains[self.ring_c + 1] = n[0] + "Br"
                    elif n[1:] == "Bz":
                        self.side_chains[self.ring_c + 1] = n[0] + "C(=O)c2ccccc2"
                    elif n[1:] == "Gc":
                        self.side_chains[self.ring_c + 1] = n[0] + "C(=O)CO"
                    elif n[1:] == "Me":
                        self.side_chains[self.ring_c + 1] = n[0] + "C"
                    elif n[1:] == "Ph":
                        self.side_chains[self.ring_c + 1] = n[0] + "c2ccccc2"
                    elif n[1:] == "Tf":
                        self.side_chains[self.ring_c + 1] = n[0] + "S(=O)(=O)C(F)(F)F"
                    elif n[1:] == "Tr":
                        self.side_chains[self.ring_c + 1] = n[0] + "C(c2ccccc2)(c3ccccc3)c4ccccc4"
                    elif n[1:] == "Ts":
                        self.side_chains[self.ring_c + 1] = n[0] + "S(=O)(=O)c2ccc(C)cc2"
                    else:
                        not_implemented_message(n)
                        full = False
                else:
                    not_implemented_message(n)
                    full = False
            elif len(n) == 4:
                if n[0].isdigit():
                    if n[1:] == "NAc":
                        self.side_chains[int(n[0])] = "NC(C)(=O)"
                    elif n[1:] == "Ala":
                        elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol()
                        self.side_chains[int(n[0])] = elem + "C(=O)[C@@H](C)N"
                    else:
                        not_implemented_message(n)
                        full = False
                    if n[1] in "ON":
                        if n[2:] == "Me":
                            self.side_chains[int(n[0])] = n[1] + "C"
                        elif n[2:] == "Ac":
                            self.side_chains[int(n[0])] = n[1] + "C(C)(=O)"
                        elif n[2:] == "Cl":
                            self.side_chains[int(n[0])] = n[1] + "Cl"
                        elif n[2:] == "Bn":
                            self.side_chains[int(n[0])] = n[1] + "c2ccccc2"
                        elif n[2:] == "Br":
                            self.side_chains[int(n[0])] = n[1] + "Br"
                        elif n[2:] == "Bz":
                            self.side_chains[int(n[0])] = n[1] + "C(=O)c2ccccc2"
                        elif n[2:] == "Gc":
                            self.side_chains[int(n[0])] = n[1] + "C(=O)CO"
                        elif n[2:] == "Ph":
                            self.side_chains[int(n[0])] = n[1] + "c2ccccc2"
                        elif n[2:] == "Tf":
                            self.side_chains[int(n[0])] = n[1] + "S(=O)(=O)C(F)(F)F"
                        elif n[2:] == "Tr":
                            self.side_chains[int(n[0])] = n[1] + "C(c2ccccc2)(c3ccccc3)c4ccccc4"
                        elif n[2:] == "Ts":
                            self.side_chains[int(n[0])] = n[1] + "S(=O)(=O)c2ccc(C)cc2"
                        else:
                            not_implemented_message(n)
                            full = False
            elif len(n) == 6:
                if n[-2] == "P":
                    self.side_chains[int(n[0])] = "OP(=O)(O)O"
                else:
                    not_implemented_message(n)
                    full = False
            elif len(n) == 7:
                if n[-3:] == "Me-":
                    self.side_chains[int(n[0])] = "OC"
                elif n[-3:] == "Ac-":
                    elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol()
                    self.side_chains[int(n[0])] = elem + "C(C)(=O)"
                else:
                    not_implemented_message(n)
                    full = False
            else:
                not_implemented_message(n)
                full = False

        self.assemble_chains()

        return self.monomer, full

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
                    if chain == "H":  # Account for desoxygenation
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
        smiles = smiles
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
