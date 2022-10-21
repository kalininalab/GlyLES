import copy
import logging
import re

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdmolops import AddHs, RemoveHs

from glyles.glycans.factory.factory_o import OpenFactory
from glyles.glycans.mono.reactor_basic import change_base_monomer
from glyles.glycans.utils import Enantiomer, ketoses2
from glyles.grammar.GlycanLexer import GlycanLexer

O, C = 0, 1

# map of all functional groups from IUPAC name to SMILES string
functional_groups = {
    "": "",
    "Aep": "NCCP([O-])(=O)[O-]",
    "Ala": "N[C@@H](C)C(=O)O",
    "Am": "C(=N)C",
    "Asp": "N[C@@H](CC(=O)O)C(=O)O",
    "Br": "Br",
    "Cer": "OC[C@H](NC=O)C(O)/C=C/" + "C" * 13,
    "Cho": "OCC[N+](C)(C)C",
    "Cl": "Cl",
    "Cm": "NC(=O)",
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
    "P": "OP(=O)(O)O",
    "PhNO2": "Oc2ccc([N+]([O-])=O)cc2",
    "Pp": "OC(=O)CC",
    "S": "S(=O)(=O)O",
    "Pro": "N2CCCC2C(=O)O",
    "Ser": "OC(=O)[C@H](N)CO",
    "Tf": "OS(=O)(=O)C(F)(F)F",
    "Thr": "N[C@@H](C(=O)O)[C@@H](O)C",
    "Ulo": "OC(c2ccccc2)(c2ccccc2(Cl))CCN(C)C",

    # COH land
    "A": "OC(=O)C",
    "Allyl": "CC=C",
    "Ac": "C(=O)C",
    "Ang": "OC(=O)C/(C)=C\C",
    "Bz": "C(=O)c2ccccc2",
    "Bn": "Cc2ccccc2",
    "cdPam": "OC(=O)CCCCCCC/C=C\CCCCCC",
    "Cet": "CCC(=O)O",
    "Cin": "OC(=O)/C=C/c2ccccc2",
    "Coum": "OC(=O)/C=C/c1ccc(O)ccc",
    "Dce": "OC(=O)CCCCCCCC=C",
    "Dhp": "OC(=O)C(O)(O)CCC",
    "Dhpa": "OC(=O)C(O)(O)CCC",
    "Etg": "OCCO",
    "Fer": "OC(=O)/C=C/c2ccc(O)c(OC)c2",
    "Fo": "OC(=O)",
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
    "Sin": "OC(=O)/C=C/c2cc(OC)c(O)c(OC)c2",
    "Suc": "OC(=O)CCC(=O)O",
    "Tig": "OC(=O)C/(C)=C/C",
    "Tr": "C(c2ccccc2)(c3ccccc3)c4ccccc4",
    "Ts": "OS(=O)(=O)c2ccc(C)cc2",
    "Vac": "OC(=O)CCCCCCCCCC=CCCCCCC",

    # Sugar rings
    "Pen": "OC2C(O)C(O)C(O)CO2",
    "Hex": "OC2C(O)C(O)C(O)C(CO)O2",
    "Hep": "OC2C(O)C(O)C(O)C(C(O)CO)O2",
    "Oct": "OC2C(O)C(O)C(O)C(C(O)C(O)CO)O2",

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
    "Phthi": "OC(=O)C(C)=CC(C)CC(C)" + "C" * 18,
    "Ner": "OC(=O)" + "C" * 13 + "/C=C\C" + "C" * 7,
}

# list of functional groups that preserve the atom it is attached to instead of replacing it
preserve_elem = [
    "Ac", "Allyl", "Am", "Bn", "Bz", "Gc", "P", "S"
]

# list of functional groups that start with an N and might be confused with a nitrogen-bridge
n_conflict = [
    "Ner", "Nno", "Non", "Nn"
]

# list of functional groups that start with an O and might be confused with an oxygen-bridge
o_conflict = [
    "Oco", "Ole", "Orn", "Oc"
]

# list of functional groups that start with a P and might be confused with a phosphate-bridge
p_conflict = [
    "Psyllic", "Prop", "Pam", "Pro", "Pyr", "Pe", "Ph", "Pr", "Pp"
]

# list of functional groups that start with a C and might be confused with a carbon-bridge
c_conflict = [
    "Cct", "Cer", "Cet", "Cho", "Cin", "Crt", "Cys", "Cl", "Cm"
]


def not_implemented_message(mod):
    """
    Print a message if some modification is parsable but not implemented yet. So, the grammar might accept some
    modifications but the reactor is yet not capable to attach that functional group to the monomer.

    Args
        mod (str): name of the modification that cannot be converted

    Returns
        Nothing
    """
    logging.warning(
        f"ModificationNotImplementedWarning: {mod} Modification not implemented. "
        f"The returned molecule will not have this modification."
    )


def extract_bridge(n):
    """
    Extract bridge and functional groups from the IUPAC description.

    Args:
        n (str): name of functional group in IUPAC format

    Returns:
        Tuple of bridge element and the functional group
    """
    # Find the lowercase part that is for sure part of the functional group's name
    i = re.search('[a-z]', n)
    if i is not None:
        i = i.start()
        if i >= 2 and n[i - 2] == "H":  # in case of something like Homoserine (HSer)
            i -= 1
        fg = n[i - 1:]
        bridge = n[:i - 1]

        # if bridge overlaps with functional group
        bridge = bridge[:-1] if functional_groups[fg][0] == bridge[-1] else bridge
    else:
        # if no lowercase in input the last char is the functional group, everything before the bridge
        fg = n[-1]
        bridge = n[:-1]

    # if bridge starts with a number, drop it
    if len(bridge) > 0 and bridge[0].isdigit():
        bridge = bridge[1:]

    # if bridge starts with carbon, drop it
    if len(bridge) > 0 and bridge[0] == "C":
        bridge = bridge[1:]

    return bridge, fg


def opposite_chirality(tag):
    """
    Invert chirality based on the given ChiralType.

    Args:
        tag (rdkit.Chem.rdchem.ChiralType): ChiralType to be inverted

    Returns:
        Returns opposite chiral tag
    """
    if tag == ChiralType.CHI_TETRAHEDRAL_CCW:
        return ChiralType.CHI_TETRAHEDRAL_CW
    elif tag == ChiralType.CHI_TETRAHEDRAL_CW:
        return ChiralType.CHI_TETRAHEDRAL_CCW
    return tag


class SMILESReaktor:
    def __init__(self, monomer):
        """
        Initialize the reaktor with the monosaccharide to attach functional groups to.

        Args:
            monomer (Monomer): monosaccharide object
        """
        self.monomer = monomer
        self.monomer.get_structure()
        self.side_chains = None
        self.ring_c = 2 if (monomer, monomer.get_lactole()) in ketoses2 else 1

    def react(self, names, types):
        """
        Main function to add functional groups to the monomer stored in this object.

        Args:
            names (List[str]): List of the names of the functional groups
            types (List[int]): List of type definitions from the recipe of a monosaccharide

        Returns:
            The monomer object and a bool flag indicating that all functional groups were recognized and attached
        """
        full = True

        self.monomer = change_base_monomer(self.monomer, names, types)

        # check for anhydro forms in the recipe
        self.check_for_anhydro(names, types)

        # store all functional groups that cannot be attached in the current round for the next round
        higher_order_groups = [], []
        while len(names) > 0:
            start_len = len(names)

            # define list of functional groups to be attached
            self.side_chains = [["", ""] for _ in range(1 + np.count_nonzero(self.monomer.x[:, 0] == 6))]

            # parse remaining modifications
            for n, t in zip(names, types):
                # if it's not a modification or already parsed, continue
                if t != GlycanLexer.MOD or \
                        n.count("L") + n.count("D") == len(n) or \
                        n in ['-', '-ol', '-onic', '-aric', '-ulosonic', '-ulosaric'] or \
                        "Anhydro" in n:
                    continue

                # if the name starts with a - char, drop this, the functional group will be independent of the previous
                if n[0] == "-":
                    n = n[1:]

                # parse making the monosaccharide an acid
                if n == "A" or n == "-uronic":
                    # if there's no ring take carbon with the highest number
                    if sum(self.monomer.x[:, 2] & 0b1) == 0:
                        c_id = int(max(self.monomer.x[self.monomer.x[:, 0] == 6, 1]))

                    # else take the last carbon in the ring
                    else:
                        c_id = int(max(self.monomer.x[(self.monomer.x[:, 0] == 6) & (self.monomer.x[:, 2] & 0b1), 1]))
                        c_id = int(np.where(self.monomer.x[:, 1] == c_id)[0])

                    # if the selected carbon has a tail raging away from the monomer, iterate all the way down
                    children = np.where((self.monomer.adjacency[c_id, :] == 1) & (self.monomer.x[:, 0] == 6) &
                                        (1 - self.monomer.x[:, 2] & 0b1))[0].tolist()
                    while len(children) != 0:
                        c_id = int(children[0])
                        children = np.where(
                            (self.monomer.adjacency[c_id, :] == 1) & (self.monomer.x[:, 0] == 6) &
                            (1 - self.monomer.x[:, 2] & 0b1) & (self.monomer.x[:, 1] > self.monomer.x[c_id, 1])
                        )[0].tolist()

                    # attach that acid group in one way or another
                    if len(self.side_chains[self.monomer.x[c_id, 1]][O]) > 0 and \
                            self.side_chains[self.monomer.x[c_id, 1]][O][-1] == "O":
                        self.side_chains[self.monomer.x[c_id, 1]][O] += "C(=O)O"
                    else:
                        self.side_chains[self.monomer.x[c_id, 1]][O] += "(=O)O"

                # add a nitrogen
                elif n == "N":
                    self.side_chains[1 if self.monomer.get_name() in ['Fru', 'Tag', 'Sor', 'Psi'] else 2][O] += "N"

                # change the enantiomer to the D or L form
                elif n == "D-":
                    self.to_enantiomer(Enantiomer.D)
                elif n == "L-":
                    self.to_enantiomer(Enantiomer.L)

                # add acid or glycoly groups to neuraminic acid
                elif n == "Ac" and self.monomer.get_name() == "Neu":
                    self.side_chains[5][O] += "NC(=O)C"
                elif n == "Gc" and self.monomer.get_name() == "Neu":
                    self.side_chains[5][O] += "NC(=O)CO"

                # if the side chain starts with a position specification
                elif n[0].isdigit():

                    # if the functional group cannot be attached at the moment, wait for next round
                    if int(n[0]) > len(self.side_chains) - 1:
                        higher_order_groups[0].append(n)
                        higher_order_groups[1].append(t)
                        continue

                    # desoxygenation of a specific position
                    if n[1:] == "d":
                        self.side_chains[int(n[0])][O] += "H"

                    # change chirality at a single chiral carbon atom
                    elif n[1:] == "e":
                        idx = int(np.where(self.monomer.x[:, 1] == int(n[0]))[0])
                        tag = opposite_chirality(self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag())
                        self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(tag)

                    # add a functional group connected with an oxygen or nitrogen
                    elif len(n) > 4 and n[1] == n[3] == "-" and n[2] in "ON":
                        elem = "" if functional_groups[n[4:-1]][0] == n[2] else n[2]
                        full &= self.set_fg(O, int(n[0]), elem, n[4:-1])

                    # connect a functional group with nitrogen, oxygen, or phosphate in between
                    elif n[1] in "NOP" and n[1:] not in n_conflict + o_conflict + p_conflict:
                        bridge, fg = extract_bridge(n)
                        if len(bridge) > 0 and functional_groups[fg][0] == bridge[-1]:
                            bridge = bridge[:-1]
                        full &= self.set_fg(O, int(n[0]), bridge, fg)

                    # add a poly carbon group
                    elif (n[1] == "C" and n[2:4].isnumeric()) or \
                            (n[1:3] in "aCiC" and n[3:5].isnumeric()) or \
                            (n[1:4] == "aiC" and n[4:6].isnumeric()):
                        full &= self.set_fg(O, int(n[0]), "", n)

                    # add a group connected directly to the C-Atom
                    elif n[1] == "C" and n[1:] not in c_conflict:
                        bridge, fg = extract_bridge(n)
                        full &= self.set_fg(C, int(n[0]), bridge, fg)

                    # otherwise just add the sidechain directly to the specified position
                    else:
                        elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol() \
                            if n[1:] in preserve_elem else ""
                        if elem == "C":
                            full &= self.set_fg(C, int(n[0]), "", n[1:])
                        else:
                            full &= self.set_fg(O, int(n[0]), elem, n[1:])

                # if the side-chain is connected with a nitrogen/oxygen/phosphate to the monomer
                elif n[0] in "NOP" and n not in n_conflict + o_conflict + p_conflict:
                    bridge, fg = extract_bridge(n)
                    if len(bridge) > 0 and functional_groups[fg][0] == bridge[-1]:
                        bridge = bridge[:-1]

                    # attach the monomer to the second carbon in the ring, it might be C2 or C3 (for 2-ketoses)
                    if fg == "Me":
                        full &= self.set_fg(O, self.ring_c, bridge, fg)
                    else:
                        full &= self.set_fg(O, self.ring_c + 1, bridge, fg)

                # functional groups might be directly connected to the carbon of the ring
                elif n[0] == "C" and n not in c_conflict:
                    if "=" in n or n[1:].isdigit():
                        self.set_fg(O, self.ring_c, "O", self.parse_poly_carbon(n))
                    else:  # add a group connected directly to the C-Atom
                        bridge, fg = extract_bridge(n)
                        full &= self.set_fg(C, self.ring_c, bridge, fg)

                # in all other cases just add the functional group
                else:
                    elem = self.monomer.structure.GetAtomWithIdx(self.monomer.find_oxygen(int(n[0]))).GetSymbol() \
                        if n[1:] in preserve_elem else ""
                    elem = "" if functional_groups[n][0] == elem else elem
                    full &= self.set_fg(O, self.ring_c, elem, n)

            # after storing all functional groups in a list, iterate over them and put them into the monomers structure
            self.assemble_chains()

            if len(higher_order_groups[0]) == start_len:
                break

            names = higher_order_groups[0]
            types = higher_order_groups[1]
            higher_order_groups = [], []

        return self.monomer, full and len(higher_order_groups[0]) == 0

    def set_fg(self, c_or_o, pos, bond_elem, name):
        """
        Add the recognized functional group to the array storing all functional groups.

        Args:
            c_or_o (int): Integer indicating if the functional group is attached to an oxygen or to the carbon directly
            pos (int): The position where to attach the group
            bond_elem (str): the element(s) to use for bridging between the functional group and the monomer
            name (str): THe IUPAC-condensed name of the functional group

        Returns:
            bool flag indicating that the functional group can be attached
        """
        # check if the group is known ...
        if name in functional_groups:
            if len(self.side_chains[pos][c_or_o]) > 0:
                if self.side_chains[pos][c_or_o][-1] == bond_elem:
                    bond_elem = ""
            if len(self.side_chains[pos][c_or_o]) > 0 and \
                    self.side_chains[pos][c_or_o][-1] == functional_groups[name][0] != "C":
                self.side_chains[pos][c_or_o] += bond_elem + functional_groups[name][1:]
            else:
                self.side_chains[pos][c_or_o] += bond_elem + functional_groups[name]
            return True

        # in case of a carbon chain
        elif (name[1] == "C" and name[2:4].isnumeric()) or \
                (name[1:3] in "aCiC" and name[3:5].isnumeric()) or \
                (name[1:4] == "aiC" and name[4:6].isnumeric()):
            self.side_chains[pos][c_or_o] += bond_elem + self.parse_poly_carbon(name)
            return True

        # otherwise print message and return accordingly
        else:
            not_implemented_message(name)
            return False

    def check_for_anhydro(self, names, types):
        """
        This method checks and inserts Anhydro-groups into the molecule.

        Args:
            names (List[str]): List of the names of the functional groups
            types (List[int]): List of type definitions from the recipe of a monosaccharide

        Returns:
            Nothing
        """
        for n, t in zip(names, types):
            if t == GlycanLexer.MOD and "Anhydro" in n:
                nums = [int(x) for x in re.findall(r'\d+', n)]
                if len(nums) != 2:
                    raise ValueError("Anhydro functional groups should have exaclty two numbers: X,Y-Anhydro-...")
                y_c = int(np.where(self.monomer.x[:, 1] == nums[1])[0])
                x_o = self.monomer.find_oxygen(nums[0])
                y_o = self.monomer.find_oxygen(nums[1])

                mol = Chem.EditableMol(self.monomer.structure)
                mol.AddBond(x_o, y_c, Chem.rdchem.BondType.SINGLE)
                mol.RemoveAtom(y_o)

                self.monomer.smiles = Chem.MolToSmiles(mol.GetMol())
                self.monomer.structure = None
                self.monomer.get_structure()

    def add_to_oxygen(self, chain, idx, placeholder):
        """
        Add a functional group to an oxygen.

        Args:
            chain (str): functional group to be attached
            idx (int): position where to bind the functional group
            placeholder (int): placeholder atom to be inserted at that position to replace it later in the SMILES

        Returns:
            Nothing
        """
        if chain == "H":  # account for desoxygenation
            self.monomer.structure.GetAtomWithIdx(idx).SetAtomicNum(1)
        else:
            self.monomer.structure.GetAtomWithIdx(idx).SetAtomicNum(placeholder)

    def add_to_carbon(self, chain, i, placeholder):
        """
        Add a functional group directly to a carbon.

        Args:
            chain (str): functional group to be attached
            i (int): position where to bind the functional group
            placeholder (int): placeholder atom to be inserted at that position to replace it later in the SMILES

        Returns:
            Nothing
        """
        # make hydrogens visible in order to replace one of them later
        tmp = AddHs(self.monomer.structure)

        # identify a hydrogen bound to the carbon atom to attach the functional group to
        h = None
        for n in tmp.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == i)[0])).GetNeighbors():
            if n.GetAtomicNum() == 1:
                h = n.GetIdx()
                break
        # if none were found, raise an error
        if h is None:
            raise ValueError(f"There is no oxygen, nitrogen, and hydrogen attached to C{i}! "
                             f"No functional group can be attached there!")

        # replace the hydrogen with the functional groups placeholder
        tmp.GetAtomWithIdx(h).SetAtomicNum(placeholder)
        self.monomer.structure = RemoveHs(tmp)
        if chain[0] == "O" and chain not in preserve_elem:
            self.side_chains[i][C] = self.side_chains[i][C][1:]

    def assemble_chains(self):
        """
        Attach the chains to the monosaccharide.

        Returns:
            Nothing
        """
        # if there is no functional group to attach, return
        if not "".join("".join(x) for x in self.side_chains):
            return

        placeholder = [
            [
                (0, ""), (31, r"\[GaH*\d*\]"), (32, r"\[GeH*\d*\]"), (33, r"\[AsH*\d*\]"), (34, r"\[SeH*\d*\]"),
                (49, r"\[InH*\d*\]"), (50, r"\[SnH*\d*\]"), (51, r"\[SbH*\d*\]"), (52, r"\[TeH*\d*\]"),
                (118, r"\[OgH*\d*\]"),
            ], [
                (0, ""), (81, r"\[TlH*\d*\]"), (82, r"\[PbH*\d*\]"), (83, r"\[BiH*\d*\]"), (84, r"\[PoH*\d*\]"),
                (113, r"\[NhH*\d*\]"), (114, r"\[FlH*\d*\]"), (115, r"\[McH*\d*\]"), (116, r"\[LvH*\d*\]"),
                (117, r"\[TsH*\d*\]"),
            ]
        ]

        # iterate over all functional groups, ...
        for i, (chain, c_chain) in enumerate(self.side_chains):
            # ... check if they exist and in case add the groups
            if chain:
                idx = self.monomer.find_oxygen(i)
                if self.monomer.structure.GetAtomWithIdx(idx).GetSymbol() != "C":
                    self.add_to_oxygen(chain, idx, placeholder[O][i][0])
                else:
                    self.add_to_carbon(chain, i, placeholder[C][i][0])
            if c_chain:
                self.add_to_carbon(c_chain, i, placeholder[C][i][0])

        # create SMILES from molecule with placeholders
        smiles = self.monomer.to_smiles(0, root_idx=100)

        ring_offset = len(self.monomer.ring_info) - 1
        # replace all placeholders by their actual functional group
        for i, (chain, c_chain) in enumerate(self.side_chains):
            if chain:
                # chain = re.sub("[\d+]", lambda x: str(int(x.group()) + ring_offset, chain))
                chain = re.sub('({})'.format(2), lambda x: str(int(x.group(1)) + ring_offset), chain)
                # smiles = re.sub(placeholder[O][i][1], "" if chain == "H" else chain, smiles)
                smiles = re.sub(placeholder[O][i][1], ("" if chain == "H" else chain).replace("\\", "\\\\"),
                                smiles).replace("\\\\", "\\")
                # smiles = smiles.replace(placeholder[O][i][1], "" if chain == "H" else chain)
            if c_chain:
                c_chain = re.sub('({})'.format(2), lambda x: str(int(x.group(1)) + ring_offset), c_chain)
                smiles = re.sub(placeholder[C][i][1], c_chain.replace("\\", "\\\\"), smiles).replace("\\\\", "\\")
                # smiles = smiles.replace(placeholder[C][i][1], chain)

        # update the monomer accordingly
        self.monomer.smiles = smiles.replace("()", "")
        self.monomer.structure = None
        self.monomer.get_structure()

    def to_enantiomer(self, form):
        """
        Convert the stored monomer into the specified enantiomer by switching all chiral tags in the rings carbon atoms.

        Args:
            form (Enantiomer): Enantiomer, specifying L or D form

        Returns:
            Nothing
        """
        # check if the enantiomeric requirement is not already fulfilled
        if self.monomer.get_isomer() == form:
            return

        # otherwise switch all chiral tags
        for idx in [atom.GetIdx() for atom in self.monomer.structure.GetAtoms()]:
            self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(
                opposite_chirality(
                    self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag()
                )
            )

    def parse_poly_carbon(self, name):
        """
        Parse long pure-carbon chains.

        Args:
            name (str): name of the long carbon chain with all the specifications of double bonds

        Returns:
            String for functional group based on the given name
        """
        # extract key information about the carbon string, such as terminal end and number of carbon atoms
        ante = name[1] == "a"
        iso = name[2 if ante else 1] == "i"
        count = int(re.findall(r'\d+', name)[1])

        # generate start and end of carbon strand and a plain carbon strand of remaining carbon atoms
        start = "OC(=O)"
        if ante and iso:
            end = "(C)CC"
            count -= 3
        elif ante:
            raise NotImplementedError()
        elif iso:
            end = "(C)C"
            count -= 2
        else:
            end = ""
        chain = "C" * (count - 1)

        # extract information on double bonds and rings in carbon strand and fill them in a mods list
        parts = re.findall(r'\{[\w|,]*}', name)
        if len(parts) != 0:
            mods = []
            for part in parts:
                pre = name[name.index(part) - 1]
                idx = part[1:-1].split(",")
                if pre == "c":
                    for i in idx:
                        mods.append((int(i), "2"))
                if pre == "=":
                    for i in idx:
                        if i[0] == "c":  # cis C=C bond
                            pos = int(i[1:])
                            mods.append((pos - 1, "/"))
                            mods.append((pos, "="))
                            mods.append((pos + 1, "\\"))
                        elif i[0] == "t":  # trans C=C bond
                            pos = int(i[1:])
                            mods.append((pos - 1, "/"))
                            mods.append((pos, "="))
                            mods.append((pos + 1, "/"))
                        else:
                            mods.append((int(i), "="))

            # check if only modifiable carbon atoms are contained in the modification list
            if (max(mods, key=lambda x: x[0])[0] - 1) > len(chain):
                raise ValueError(f"{name} defines modifications for carbon atoms that are not modifiable.")

            # sort the modifications in descending order based on their indices
            for pos, mod in sorted(mods, key=lambda x: x[0], reverse=True):
                pos = pos - 1  # as C1 is covered up in start
                chain = chain[:pos] + mod + chain[pos:]

        # glue everything together and return it
        return (start + chain + end).replace("//", "/")
