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
    "Am": "C(C)(=N)",
    "Cl": "Cl",
    "Bn": "OCc2ccccc2",
    "Br": "Br",
    "Bz": "C(=O)c2ccccc2",
    "Fo": "OC=O",
    "Gc": "C(=O)CO",
    "Ph": "c2ccccc2",
    "Pp": "OC(=O)CC",
    "Tf": "OS(=O)(=O)C(F)(F)F",
    "Tr": "C(c2ccccc2)(c3ccccc3)c4ccccc4",
    "Ts": "OS(=O)(=O)c2ccc(C)cc2",
    "Ang": "OC(=O)C/(C)=C\C",
    "Asp": "N[C@@H](CC(=O)O)C(=O)O",
    "Beh": "OC(=O)" + "C" * 21,
    "Cho": "OCC[N+](C)(C)C",
    "Cys": "N[C@@H](CS)C(=O)O",
    "Gly": "OCC(O)CO",
    "Lys": "NCCCC[C@H](N)C(=O)O",
    "Mal": "C(C(=O)O)CC(=O)O",
    "Ole": "OC(=O)CCCCCCC/C=C\CCCCCCCC",
    "Pro": "N1CCCC1C(=O)O",
    "Tig": "OC(=O)C/(C)=C/C",
    "Thr": "N[C@H](C(=O)O)[C@@H](O)C",
    "Dce": "OC(=O)CCCCCCCC=C",
    "3oxoMyr": "OC(=O)CC(=O)CCCCCCCCCCC",
    "dPam": "OC(=O)CCCCCCCC=CCCCCCC",
    "cdPam": "OC(=O)CCCCCCC/C=C\CCCCCC",
    "tdPam": "OC(=O)CCCCCCC/C=C/CCCCCC",
    "Vac": "OC(=O)CCCCCCCCCC=CCCCCCC",
    "Lin": "OC(=O)CCCCCCC/C=C\C/C=C\CCCCC",
    "cVac": "OC(=O)CCCCCCCCC/C=C\CCCCCC",
    "17HOLin": "OC(=O)CCCCCCC/C=C\C/C=C\CCCC(O)C",
    "aLnn": "OC(=O)CCCCCCCC=CCC=CCC=CCC",
    "gLnn": "OC(=O)CCCCC=CCC=CCC=CCCCCC",
    "eSte": "OC(=O)CCCCCCCC=CC=CC=CCCCC",
    "d2Ach": "OC(=O)CCCCCCC=CCC=CCCCCCCCC",
    "d3Ach": "OC(=O)CCCC=CCC=CCC=CCCCCCCCC",
    "d4Ach": "OC(=O)CCCC=CCC=CCC=CCC=CCCCCC",
    "NFo": "NC=O",
    "Phyt": "OCCC(C)CCCC(C)CCCC(C)CCCC(C)C",
    "Allyl": "CC=C",
    "PhNO2": "Oc1ccc([N+]([O-])=O)cc1",
    "tBu": "OC(C)(C)C",
    "But": "OC(=O)" + "C" * 2,
    "Vl": "OC(=O)" + "C" * 4,
    "Hxo": "OC(=O)" + "C" * 5,
    "Hpo": "OC(=O)" + "C" * 6,
    "Oco": "OC(=O)" + "C" * 7,
    "Nno": "OC(=O)" + "C" * 8,
    "Non": "OC(=O)" + "C" * 8,
    "Dco": "OC(=O)" + "C" * 9,
    "Udo": "OC(=O)" + "C" * 10,
    "Lau": "OC(=O)" + "C" * 11,
    "Myr": "OC(=O)" + "C" * 13,
    "Pam": "OC(=O)" + "C" * 15,
    "Mar": "OC(=O)" + "C" * 16,
    "Ste": "OC(=O)" + "C" * 17,
    "Ach": "OC(=O)" + "C" * 19,
    "Lig": "OC(=O)" + "C" * 23,
    "Crt": "OC(=O)" + "C" * 25,
    "Ccr": "OC(=O)" + "C" * 26,
    "Mon": "OC(=O)" + "C" * 27,
    "Mel": "OC(=O)" + "C" * 28,
    "Lacceroic": "OC(=O)" + "C" * 31,
    "Psyllic": "OC(=O)" + "C" * 32,
    "Geddic": "OC(=O)" + "C" * 33,
    "Ceroplastic": "OC(=O)" + "C" * 35,
    "Me": "O" + "C" * 1,
    "Et": "O" + "C" * 2,
    "Pr": "O" + "C" * 3,
    "Bu": "O" + "C" * 4,
    "Pe": "O" + "C" * 5,
    "Hx": "O" + "C" * 6,
    "Hp": "O" + "C" * 7,
    "Oc": "O" + "C" * 8,
    "Nn": "O" + "C" * 9,
    "Dec": "O" + "C" * 10,
    "Und": "O" + "C" * 11,
    "Dod": "O" + "C" * 12,
    "Fer": "OC(=O)/C=C/c1ccc(O)c(OC)c1",
    "Etg": "OCCO",
    "Cin": "OC(=O)/C=C/c1ccccc1",
    "Cm": "NC=O",
    "Etn": "OCCN",
    "EtN": "OCCN",
    "Glu": "OC(=O)CC[C@H](N)C(=O)O",
    "Hse": "OCCC(N)C(=O)O",
    "Orn": "NCCC[C@H](N)C(=O)O",
    "Pyr": "OC(=O)C(=O)C",
    "Ser": "OC[C@H](N)C(=O)O",
    "Sin": "OC(=O)/C=C/c1cc(OC)c(O)c(OC)c1",
}

preserve_elem = [
    "Ac", "Allyl", "Am", "Bz", "Gc", "P", "S"
]


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
                elif n[0] == "C" and ("=" in n or n[1:].isdigit()):
                    self.set_fg(self.ring_c, "O", self.parse_poly_carbon(n))
                else:
                    elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol() \
                        if n[1:] in preserve_elem else ""
                    full &= self.set_fg(int(n[0]), elem, n[1:])
            elif n[0] in "NO" and n not in ["Ole"]:
                if n[1:] == "Me":
                    full &= self.set_fg(self.ring_c, n[0], n[1:])
                else:
                    full &= self.set_fg(self.ring_c + 1, n[0], n[1:])
            elif n[0] == "C" and ("=" in n or n[1:].isdigit()):
                self.set_fg(self.ring_c, "O", self.parse_poly_carbon(n))
            else:
                elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol() \
                    if n[1:] in preserve_elem else ""
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

    def parse_poly_carbon(self, name):
        pass
