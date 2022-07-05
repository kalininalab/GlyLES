import copy
import logging

import numpy as np
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdmolops import AddHs, RemoveHs

from glyles.glycans.factory.factory_o import OpenFactory
from glyles.glycans.utils import Enantiomer, ketoses2
from glyles.grammar.GlycanLexer import GlycanLexer

O, C = 0, 1

functional_groups = {
    "": "",
    "Aep": "NCCP([O-])([O-])=O",
    "Ala": "N[C@@H](C)C(=O)O",
    "Am": "C(C)(=N)",
    "Asp": "N[C@@H](CC(=O)O)C(=O)O",
    "Br": "Br",
    "Cho": "OCC[N+](C)(C)C",
    "Cl": "Cl",
    "Cm": "NC=O",
    "Cys": "N[C@@H](CS)C(=O)O",
    "Etn": "OCCN",
    "EtN": "OCCN",
    "F": "F",
    "Glu": "OC(=O)CC[C@H](N)C(=O)O",
    "Hse": "OCCC(N)C(=O)O",
    "HSer": "OC(=O)[C@H](N)CCO",
    "I": "I",
    "Lys": "NCCCC[C@H](N)C(=O)O",
    "N": "N",
    "NFo": "NC=O",
    "Orn": "NCCC[C@H](N)C(=O)O",
    "P": "P(=O)(O)(O)",
    "PhNO2": "Oc1ccc([N+]([O-])=O)cc1",
    "Pp": "OC(=O)CC",
    "S": "S(=O)(=O)(O)",
    "Pro": "N1CCCC1C(=O)O",
    "Ser": "OC(=O)[C@H](N)CO",
    "Tf": "OS(=O)(=O)C(F)(F)F",
    "Thr": "N[C@H](C(=O)O)[C@@H](O)C",
    "Ulo": "OC(c1ccccc1)(c1ccccc1(Cl))CCN(C)C",

    # COH land
    "Allyl": "CC=C",
    "Ac": "C(C)(=O)",
    "Ang": "OC(=O)C/(C)=C\C",
    "Bz": "C(=O)c2ccccc2",
    "Bn": "OCc2ccccc2",
    "cdPam": "OC(=O)CCCCCCC/C=C\CCCCCC",
    "Cet": "CCC(=O)O",
    "Cin": "OC(=O)/C=C/c1ccccc1",
    "Dce": "OC(=O)CCCCCCCC=C",
    "Dhp": "OC(=O)C(O)(O)CCC",
    "Dhpa": "OC(=O)C(O)(O)CCC",
    "Etg": "OCCO",
    "Fer": "OC(=O)/C=C/c1ccc(O)c(OC)c1",
    "Fo": "OC=O",
    "Gc": "C(=O)CO",
    "Gly": "OCC(O)CO",
    "Gro": "OCC(O)CO",
    "He": "C(O)C",
    "Lac": "OC(=O)C(O)C",
    "Lin": "OC(=O)CCCCCCC/C=C\C/C=C\CCCCC",
    "Mal": "O[C@H](C(=O)O)CC(=O)O",
    "Ole": "OC(=O)CCCCCCC/C=C\CCCCCCCC",
    "Ph": "c2ccccc2",
    "Phyt": "OCCC(C)CCCC(C)CCCC(C)CCCC(C)C",
    "Pyr": "OC(=O)C(=O)C",
    "Sin": "OC(=O)/C=C/c1cc(OC)c(O)c(OC)c1",
    "Suc": "OC(=O)CCC(=O)O",
    "Tig": "OC(=O)C/(C)=C/C",
    "Tr": "C(c2ccccc2)(c3ccccc3)c4ccccc4",
    "Ts": "OS(=O)(=O)c2ccc(C)cc2",
    "Vac": "OC(=O)CCCCCCCCCC=CCCCCCC",

    # Sugar rings
    "Pen": "OC1C(O)C(O)C(O)CO1",
    "Hex": "OC1C(O)C(O)C(O)C(CO)O1",
    "Hep": "OC1C(O)C(O)C(O)C(C(O)CO)O1",
    "Oct": "OC1C(O)C(O)C(O)C(C(O)C(O)CO)O1",

    # some special cases
    "3oxoMyr": "OC(=O)CC(=O)CCCCCCCCCCC",
    "17HOLin": "OC(=O)CCCCCCC/C=C\C/C=C\CCCC(O)C",
    "aLnn": "OC(=O)CCCCCCCC=CCC=CCC=CCC",
    "cVac": "OC(=O)CCCCCCCCC/C=C\CCCCCC",
    "d2Ach": "OC(=O)CCCCCCC=CCC=CCCCCCCCC",
    "d3Ach": "OC(=O)CCCC=CCC=CCC=CCCCCCCCC",
    "d4Ach": "OC(=O)CCCC=CCC=CCC=CCC=CCCCCC",
    "dPam": "OC(=O)CCCCCCCC=CCCCCCC",
    "eSte": "OC(=O)CCCCCCCC=CC=CC=CCCCC",
    "gLnn": "OC(=O)CCCCC=CCC=CCC=CCCCCC",
    "tdPam": "OC(=O)CCCCCCC/C=C/CCCCCC",
    "tBu": "OC(C)(C)C",

    # Methanol, Ethanol, Propanol, ...
    "Me": "O" + "C" * 1,
    "Et": "O" + "C" * 2,
    "Pr": "O" + "C" * 3,
    "Prop": "O" + "C" * 3,
    "Bu": "O" + "C" * 4,
    "Pe": "O" + "C" * 5,
    "Hx": "O" + "C" * 6,
    "Hp": "O" + "C" * 7,
    "Oc": "O" + "C" * 8,
    "Nn": "O" + "C" * 9,
    "Dec": "O" + "C" * 10,
    "Und": "O" + "C" * 11,
    "Dod": "O" + "C" * 12,
    # ... and their acids
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
    "Hico": "OC(=O)" + "C" * 20,
    "Beh": "OC(=O)" + "C" * 21,
    "Lig": "OC(=O)" + "C" * 23,
    "Crt": "OC(=O)" + "C" * 25,
    "Ccr": "OC(=O)" + "C" * 26,
    "Mon": "OC(=O)" + "C" * 27,
    "Mel": "OC(=O)" + "C" * 28,
    "Lacceroic": "OC(=O)" + "C" * 31,
    "Psyllic": "OC(=O)" + "C" * 32,
    "Geddic": "OC(=O)" + "C" * 33,
    "Ceroplastic": "OC(=O)" + "C" * 35,
}

preserve_elem = [
    "Ac", "Allyl", "Am", "Bz", "Gc", "P", "S"
]

n_conflict = [
    "Nno", "Non", "Nn"
]

o_conflict = [
    "Oco", "Ole", "Orn", "Oc"
]

c_conflict = [
    "Cct", "Cer", "Cet", "Cho", "Cin", "Crt", "Cys", "Cl", "Cm"
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

        self.check_for_open_form(names, types)

        self.side_chains = [["", ""] for _ in range(1 + np.count_nonzero(self.monomer.x[:, 0] == 6))]

        # parse remaining modifications
        for n, t in zip(names, types):
            if t != GlycanLexer.MOD or n.count("L") + n.count("D") == len(n) or n in ['-', '-ol', '-onic']:
                continue
            if n == "A":
                self.side_chains[int(max(self.monomer.x[self.monomer.x[:, 0] == 6, 1]))][O] = "(=O)O"
            elif n == "N":
                self.side_chains[1 if self.monomer.get_name() in ['Fru', 'Tag', 'Sor', 'Psi'] else 2][O] = "N"
            elif n == "D-":
                self.to_enantiomer(Enantiomer.D)
            elif n == "L-":
                self.to_enantiomer(Enantiomer.L)
            elif n == "Ac" and self.monomer.get_name() == "Neu":
                self.side_chains[5][O] = "NC(C)(=O)"
            elif n == "Gc" and self.monomer.get_name() == "Neu":
                self.side_chains[5][O] = "NC(=O)CO"
            elif n[0].isdigit():
                if n[1:] == "d":  # desoxygenation of a specific position
                    self.side_chains[int(n[0])][O] = "H"
                elif n[1:] == "e":  # change chirality at a single chiral carbon atom
                    idx = int(np.where(self.monomer.x[:, 1] == int(n[0]))[0])
                    tag = opposite_chirality(self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag())
                    self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(tag)
                elif len(n) > 4 and n[1] == n[3] == "-" and n[2] in "ON":  # add a functional group connected with an oxygen
                    elem = "" if functional_groups[n[4:-1]][0] == n[2] else n[2]
                    full &= self.set_fg(O, int(n[0]), elem, n[4:-1])
                elif n[1] in "ON" and n[1:] not in n_conflict + o_conflict:  # connect a functional group with N or O in between
                    elem = "" if n[2:] != "" and functional_groups[n[2:]][0] == n[1] else n[1]
                    full &= self.set_fg(O, int(n[0]), elem, n[2:])
                elif n[1] == "C" and n[1:] not in c_conflict:  # add a group connected directly to the C-Atom
                    full &= self.set_fg(C, int(n[0]), "", n[2:])
                else:
                    elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol() \
                        if n[1:] in preserve_elem else ""
                    elem = "" if elem == "C" else elem
                    full &= self.set_fg(O, int(n[0]), elem, n[1:])
            elif n[0] in "NO" and n not in n_conflict + o_conflict:
                elem = "" if functional_groups[n[1:]][0] == n[0] else n[0]
                if n[1:] == "Me":
                    full &= self.set_fg(O, self.ring_c, elem, n[1:])
                else:
                    full &= self.set_fg(O, self.ring_c + 1, elem, n[1:])
            elif n[0] == "C" and n not in c_conflict:
                if "=" in n or n[1:].isdigit():
                    self.set_fg(O, self.ring_c, "O", self.parse_poly_carbon(n))
                else:  # add a group connected directly to the C-Atom
                    full &= self.set_fg(C, self.ring_c, "", n[1:])
            else:
                elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol() \
                    if n[1:] in preserve_elem else ""
                elem = "" if functional_groups[n][0] == elem else elem
                full &= self.set_fg(O, self.ring_c, elem, n)

        self.assemble_chains()

        return self.monomer, full

    def set_fg(self, c_or_o, pos, bond_elem, name):
        if name in functional_groups:
            self.side_chains[pos][c_or_o] = bond_elem + functional_groups[name]
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
        if sac_index != 0 and names[sac_index - 1].count("L") + \
                names[sac_index - 1].count("D") == len(names[sac_index - 1]):
            orientations = names[sac_index - 1]
        else:
            orientations = ""

        c_count = np.count_nonzero(self.monomer.x[:, 0] == 6)
        if names[len_index] == "Pen":
            count = 5
        elif names[len_index] == "Hex":
            count = 6
        elif names[len_index] == "Hep":
            count = 7
        elif names[len_index] == "Oct":
            count = 8
        else:
            count = c_count
        extension = "[C?H](O)" * (count - c_count) + "CO"

        for c in orientations[:-1]:
            if c == "L":
                extension = extension.replace("[C?H]", "[C@@H]", 1)
            if c == "D":
                extension = extension.replace("[C?H]", "[C@H]", 1)
        extension = extension.replace("[C?H]", "C")

        # find oxygen of c6 to delete it and find c6 to replace it
        c_id = int(np.where(self.monomer.x[:, 1] == int(max(self.monomer.x[self.monomer.x[:, 0] == 6, 1])))[0])
        # if the carbon to be extended is part of the ring, put the extension into a side_chain
        if self.monomer.x[c_id, 2] == 1:
            extension = extension.replace(")", ")(", 1) + ")"
        # if the original molecule has no oxygen at its last carbon, remove the first oxygen from the extension
        if self.monomer.x[((self.monomer.adjacency[c_id] - self.monomer.x[:, 2]) > 0) & (self.monomer.x[:, 0] == 8), 0].size == 0:
            extension = extension.replace("(O)", "", 1)
        self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(c_count)).SetAtomicNum(50)
        self.monomer.structure.GetAtomWithIdx(c_id).SetAtomicNum(32)

        # generate SMILES from rings oxygen and replace c6's oxygen and c6 by the extension
        smiles = self.monomer.to_smiles(0, root_idx=10)
        smiles = smiles.replace("[SnH]", "").replace("[GeH2]", extension).replace("[GeH]", extension).replace("()", "")

        # update the monomer
        self.monomer.smiles = smiles
        self.monomer.structure = None
        self.monomer.get_structure()

    def check_for_open_form(self, names, types):
        if not (("-ol" in names) ^ ("-onic" in names)):
            return

        sac_index = types.index(GlycanLexer.SAC)
        params = copy.copy(OpenFactory()[names[sac_index] + "-ol"])
        if "-onic" in names:
            params["smiles"] = params["smiles"].replace("C", "C(=O)", 1)

        self.monomer.smiles = params["smiles"]
        self.monomer.c1_find = params["c1_find"]
        self.monomer.structure = None
        self.monomer.get_structure()

    def add_to_oxygen(self, chain, idx, placeholder):
        if chain == "H":  # account for desoxygenation
            self.monomer.structure.GetAtomWithIdx(idx).SetAtomicNum(1)
        else:
            self.monomer.structure.GetAtomWithIdx(idx).SetAtomicNum(placeholder)

    def add_to_carbon(self, chain, i, placeholder):
        tmp = AddHs(self.monomer.structure)
        h = None
        for n in tmp.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == i)[0])).GetNeighbors():
            if n.GetAtomicNum() == 1:
                h = n.GetIdx()
                break
        if h is None:
            raise ValueError(f"There is no oxygen, nitrogen, and hydrogen attached to C{i}! "
                             f"No functional group can be attached there!")
        tmp.GetAtomWithIdx(h).SetAtomicNum(placeholder)
        self.monomer.structure = RemoveHs(tmp)
        if chain[0] == "O" and chain not in preserve_elem:
            self.side_chains[i][C] = self.side_chains[i][C][1:]

    def assemble_chains(self):
        placeholder = [
            (31, "[GaH2]"), (32, "[GeH3]"),
            (49, "[InH2]"), (50, "[SnH]"), (51, "[SbH2]"), (52, "[TeH]"),
            (81, "[TlH2]"), (82, "[PbH]"), (83, "[BiH2]"), (84, "[PoH]"),
        ]
        for i, (chain, c_chain) in enumerate(self.side_chains):
            if chain:
                idx = self.monomer.find_oxygen(i)
                if self.monomer.structure.GetAtomWithIdx(idx).GetSymbol() != "C":
                    self.add_to_oxygen(chain, idx, placeholder[i][0])
                else:
                    self.add_to_carbon(chain, i, placeholder[i][0])
        smiles = self.monomer.to_smiles(0, root_idx=10)
        for i, (chain, c_chain) in enumerate(self.side_chains):
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