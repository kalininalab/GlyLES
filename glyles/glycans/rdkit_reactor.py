import numpy as np
from rdkit import Chem
from rdkit.Chem.rdchem import Atom, BondType, EditableMol, ChiralType
import logging

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


def ringC(monomer):
    return 2 if (monomer, monomer.get_lactole()) in ketoses2 else 1


class Reactor:
    """
    Class with access to protected classes fields managing the modifications of a root monomer.
    #https://pubchem.ncbi.nlm.nih.gov/compound/N-Glycolyl-Neuraminic-acid
    """

    def __init__(self, monomer):
        """
        Initialize the reactor with the monomer to modify.

        Args:
            monomer (RDKitMonomer): monomer to be modified
        """
        self.monomer = monomer
        if self.monomer.structure is None:
            self.monomer.get_structure()

    def react(self, names, types):
        """
        Manage the parsed modifications and apply them in turn.
        Untested modifications are marked with >TBT< - to be tested in the code

        Args:
            names (List[str]): name (string representation) of the modification
            types (List[int]): Type of the parsed stings based on GlycanLexer.TYPE

        Returns:
            Modified monomer
        """
        full = True
        for name, t in zip(names, types):
            if t != GlycanLexer.MOD:
                continue
            if len(name) == 1:
                if name[0] == "N":  # put a nitrogen at position 2
                    self.set_nitrogen()
                if name[0] == "A":  # put an acid group at position 6
                    self.make_acid()
            elif len(name) == 2:
                if name[0].isdigit() or name[0] == "O":
                    if name[1] == "d":  # ?d - desoxygenate some positions in monomers  (TBT)
                        self.desoxygenate(int(name[0]))
                    # having a fluor atom instead of a hydrogen opposite to an oxygen at a certain position (TBI)
                    elif name[1] == "F":
                        not_implemented_message(name)
                        full = False
                    elif name[1] == "N":  # add a nitrogen atom to a certain position (TBT with O)
                        self.set_nitrogen(position=(ringC(self.monomer) if name[0] == "O" else int(name[0])))
                    elif name[1] == "S":  # add a sulfur atom to a certain position (TBT with O)
                        self.add_sulfur(position=(ringC(self.monomer) if name[0] == "O" else int(name[0])))
                    elif name[1] == "P":  # add a phosphate atom to a certain position (TBT with O)
                        self.add_phosphate(position=(ringC(self.monomer) if name[0] == "O" else int(name[0])))
                else:
                    if name == "D-":  # have the monomer in D form (regarding the enantiomerism)
                        self.to_enantiomer(Enantiomer.D)
                    elif name == "L-":  # have the monomer is L form (regarding the enantiomerism)
                        self.to_enantiomer(Enantiomer.L)
            elif len(name) == 3:
                if name[0].isdigit() or name[0] == "O":
                    # add a methyl group to a certain position (or to an oxygen at position 2)
                    if name.endswith("Me"):
                        self.add_methyl(position=(ringC(self.monomer) if name[0] == "O" else int(name[0])))
                    # add an acid group to a certain position (or to an oxygen at position 2)
                    elif name.endswith("Ac"):
                        self.add_acid(position=(ringC(self.monomer) if name[0] == "O" else int(name[0])))
                    # add a benzoyl group to a certain position (or to an oxygen at position 2)
                    elif name.endswith("Bn"):
                        not_implemented_message(name)
                        full = False
                        # self.add_benzyl(position=("O" if name[0] == "O" else int(name[0])))
                    elif name.endswith("Bz"):
                        not_implemented_message(name)
                        full = False
                        # self.add_benzoyl(position=("O" if name[0] == "O" else int(name[0])))
                    # add a glycolyl group to a certain position (or to an oxygen at position 2)
                    elif name.endswith("Gc"):
                        not_implemented_message(name)
                        full = False
                        # self.add_glycolyl(position=("O" if name[0] == "O" else int(name[0])))
                    elif name.endswith("Ph"):
                        self.add_phenyl(position=(ringC(self.monomer) if name[0] == "O" else int(name[0])))
                    elif name.endswith("Tf"):
                        not_implemented_message(name)
                        full = False
                        # self.add_triflyl(position=("O" if name[0] == "O" else int(name[0])))
                    elif name.endswith("Tr"):
                        not_implemented_message(name)
                        full = False
                        # self.add_trityl(position=("O" if name[0] == "O" else int(name[0])))
                    elif name.endswith("Ts"):
                        not_implemented_message(name)
                        full = False
                        # self.add_tosyl(position=("O" if name[0] == "O" else int(name[0])))
                elif name == "NAc":  # add a nitrogen and attach an acid group to that at position 2
                    self.add_acid(pos=self.set_nitrogen())
                elif name == "NBz":  # add a nitrogen and attach a benzoyl group to that at position ?(2)
                    not_implemented_message(name)
                    full = False
                    # self.add_benzoyl(pos=self.set_nitrogen())
            elif len(name) == 7:
                if name.endswith("Me-"):  # add a methyl group to a certain position
                    self.add_methyl(position=int(name[0]))
                if name.endswith("Ac-"):  # add an acid group to a certain position
                    self.add_acid(position=int(name[0]))
                if name.endswith("Bn-"):
                    not_implemented_message(name)
                    full = False
                    # self.add_benzyl(position=int(name[0]))
                if name.endswith("Bz-"):
                    not_implemented_message(name)
                    full = False
                    # self.add_benzoyl(position=int(name[0]))
                if name.endswith("Et-"):
                    not_implemented_message(name)
                    full = False
                    # self.add_ethyl(position=int(name[0]))
            elif len(name) == 12:  # ?,?-Anhydro- ??
                not_implemented_message(name)
                full = False
                # self.add_anhydro(int(name[0]), int(name[2]))
            else:
                not_implemented_message(name)
                full = False

        return self.monomer, full

    def add_sulfur(self, position):
        """
        Add a SO3- group at the oxygen of the specified position.
        Example: Gal -> Gal3S

        Args:
            position (int): id of the carbon where to add the SO3- group to the bound oxygen

        Returns:
            Nothing
        """
        pos = self.monomer.find_oxygen(position)
        emol = EditableMol(self.monomer.structure)

        s_id = EditableMol.AddAtom(emol, Atom(16))
        o1_id = EditableMol.AddAtom(emol, Atom(8))
        o2_id = EditableMol.AddAtom(emol, Atom(8))
        o3_id = EditableMol.AddAtom(emol, Atom(8))
        EditableMol.AddBond(emol, s_id, o1_id, order=BondType.DOUBLE)
        EditableMol.AddBond(emol, s_id, o2_id, order=BondType.SINGLE)
        EditableMol.AddBond(emol, s_id, o3_id, order=BondType.DOUBLE)
        EditableMol.AddBond(emol, pos, s_id, order=BondType.SINGLE)

        self.monomer.structure = emol.GetMol()

        new_x, new_adj = self._extend_matrices(4)
        new_x[s_id:, 0] = [16, 8, 8, 8]
        self.monomer.x = new_x

        new_adj[s_id, o1_id] = 1
        new_adj[s_id, o2_id] = 1
        new_adj[s_id, o3_id] = 1
        new_adj[s_id, pos] = 1
        new_adj[o1_id, s_id] = 1
        new_adj[o2_id, s_id] = 1
        new_adj[o3_id, s_id] = 1
        new_adj[pos, s_id] = 1

        self.monomer.adjacency = new_adj

    def add_phosphate(self, position):
        """
        Add a PO3 group at the oxygen of the specified position.
        Example: Gal -> Gal3P

        Args:
            position (int): id of the carbon where to add the PO3 group to the bound oxygen

        Returns:
            Nothing
        """

        pos = self.monomer.find_oxygen(position)
        emol = EditableMol(self.monomer.structure)

        s_id = EditableMol.AddAtom(emol, Atom(15))
        o1_id = EditableMol.AddAtom(emol, Atom(8))
        o2_id = EditableMol.AddAtom(emol, Atom(8))
        o3_id = EditableMol.AddAtom(emol, Atom(8))
        EditableMol.AddBond(emol, s_id, o1_id, order=BondType.SINGLE)
        EditableMol.AddBond(emol, s_id, o2_id, order=BondType.SINGLE)
        EditableMol.AddBond(emol, s_id, o3_id, order=BondType.DOUBLE)
        EditableMol.AddBond(emol, pos, s_id, order=BondType.SINGLE)

        self.monomer.structure = emol.GetMol()

        new_x, new_adj = self._extend_matrices(4)
        new_x[s_id:, 0] = [15, 8, 8, 8]
        self.monomer.x = new_x

        new_adj[s_id, o1_id] = 1
        new_adj[s_id, o2_id] = 1
        new_adj[s_id, o3_id] = 1
        new_adj[s_id, pos] = 1
        new_adj[o1_id, s_id] = 1
        new_adj[o2_id, s_id] = 1
        new_adj[o3_id, s_id] = 1
        new_adj[pos, s_id] = 1

        self.monomer.adjacency = new_adj

    def add_acid(self, position=None, pos=None):
        """
        Add an acid group to a specific position. Here the position can be provided either as the C-index
        (position) or as the rdkit id of the atom where to append the acid group (implemented for NAc). Exactly one
        of both arguments must be provided.
        Example: GalN -> GalNAc or Gal -> Gal5Ac

        Args:
            position (int): index of the c-atom where to append the acid group
            pos (int): rdkit id of the atom where to append the acid group

        Returns:
            Nothing
        """
        if (position is None) == (pos is None):
            raise ValueError()

        if position is not None:
            pos = self.monomer.find_oxygen(position)

        emol = EditableMol(self.monomer.structure)

        c1_id = EditableMol.AddAtom(emol, Atom(6))
        c2_id = EditableMol.AddAtom(emol, Atom(6))
        o1_id = EditableMol.AddAtom(emol, Atom(8))
        EditableMol.AddBond(emol, c1_id, pos, order=BondType.SINGLE)
        EditableMol.AddBond(emol, c2_id, c1_id, order=BondType.SINGLE)
        EditableMol.AddBond(emol, o1_id, c1_id, order=BondType.DOUBLE)

        self.monomer.structure = emol.GetMol()

        new_x, new_adj = self._extend_matrices(3)
        new_x[c1_id:, 0] = [6, 6, 8]
        self.monomer.x = new_x

        new_adj[c1_id, c2_id] = 1
        new_adj[c1_id, o1_id] = 1
        new_adj[c1_id, pos] = 1
        new_adj[c2_id, c1_id] = 1
        new_adj[o1_id, c1_id] = 1
        new_adj[pos, c1_id] = 1

        self.monomer.adjacency = new_adj

    def add_methyl(self, position=None, pos=None):
        """
        Append a methyl group to the monomer at the specified position.

        Args:
            position (int): index of the c-atom where to append the acid group
            pos (int): rdkit id of the atom where to append the acid group

        Returns:
            Nothing
        """

        if (position is None) == (pos is None):
            raise ValueError()

        if position is not None:
            pos = self.monomer.find_oxygen(position)

        emol = EditableMol(self.monomer.structure)

        c_id = EditableMol.AddAtom(emol, Atom(6))
        EditableMol.AddBond(emol, c_id, pos, order=BondType.SINGLE)

        self.monomer.structure = emol.GetMol()

        new_x, new_adj = self._extend_matrices(1)
        new_x[c_id, 0] = 6
        self.monomer.x = new_x

        new_adj[c_id, pos] = 1
        new_adj[pos, c_id] = 1

        self.monomer.adjacency = new_adj

    def add_benzoyl(self, position=None, pos=None):
        """
        Add a benzoyl group to the monomer at the specified position

        Args:
            position (int): index of the c-atom where to append the benzoyl group
            pos (int): rdkit id of the atom where to append the acid group

        Returns:
            Nothing
        """
        # TODO: Add the benzoyl group
        pass

    def add_glycolyl(self, position=2):
        """
        Add a glycolyl group to the monomer at the specified position

        Args:
            position (int): Position where to add a benzoyl group

        Returns:
            Nothing
        """
        # TODO: Add the glycolyl group
        pass

    def add_phenyl(self, position=None, pos=None):
        """
        Add a phenyl group to the monomer at the specified position

        Args:
            position (int): index of the c-atom where to append the phenyl group
            pos (int): rdkit id of the atom where to append the phenyl group

        Returns:
            Nothing
        """
        print()
        if (position is None) == (pos is None):
            raise ValueError()

        if position is not None:
            print(Chem.MolToSmiles(self.monomer.structure))
            chirality = self.desoxygenate(position)
            print("Chi1:", chirality)
            pos = int(np.where(self.monomer.x[:, 1] == position)[0])
            pos2 = int(np.where(self.monomer.x[:, 1] == (position + 1))[0])
            print("Chi2:", self.monomer.structure.GetAtomWithIdx(pos2).GetChiralTag())
            print(Chem.MolToSmiles(self.monomer.structure))
        else:
            chirality = None

        emol = EditableMol(self.monomer.structure)

        c1_id = EditableMol.AddAtom(emol, Atom(6))
        c2_id = EditableMol.AddAtom(emol, Atom(6))
        c3_id = EditableMol.AddAtom(emol, Atom(6))
        c4_id = EditableMol.AddAtom(emol, Atom(6))
        c5_id = EditableMol.AddAtom(emol, Atom(6))
        c6_id = EditableMol.AddAtom(emol, Atom(6))
        EditableMol.AddBond(emol, pos, c1_id, order=BondType.SINGLE)
        EditableMol.AddBond(emol, c1_id, c2_id, order=BondType.AROMATIC)
        EditableMol.AddBond(emol, c2_id, c3_id, order=BondType.AROMATIC)
        EditableMol.AddBond(emol, c3_id, c4_id, order=BondType.AROMATIC)
        EditableMol.AddBond(emol, c4_id, c5_id, order=BondType.AROMATIC)
        EditableMol.AddBond(emol, c5_id, c6_id, order=BondType.AROMATIC)
        EditableMol.AddBond(emol, c6_id, c1_id, order=BondType.AROMATIC)

        self.monomer.structure = emol.GetMol()
        print(Chem.MolToSmiles(self.monomer.structure))
        print("Chi2:", self.monomer.structure.GetAtomWithIdx(pos2).GetChiralTag())
        print("Chi1:", self.monomer.structure.GetAtomWithIdx(pos).GetChiralTag())
        print("P1-N:", [atom.GetIdx() for atom in self.monomer.structure.GetAtomWithIdx(pos).GetNeighbors()])
        self.monomer.structure.GetAtomWithIdx(pos).SetChiralTag(chirality)
        print("Chi1:", self.monomer.structure.GetAtomWithIdx(pos).GetChiralTag())
        print(Chem.MolToSmiles(self.monomer.structure))

        new_x, new_adj = self._extend_matrices(6)
        new_x[c1_id:, 0] = [6, 6, 6, 6, 6, 6]
        new_x[c1_id:, 2] = [2, 2, 2, 2, 2, 2]
        self.monomer.x = new_x

        for c1 in [c1_id, c2_id, c3_id, c4_id, c5_id, c6_id]:
            for c2 in [c1_id, c2_id, c3_id, c4_id, c5_id, c6_id]:
                if c1 != c2:
                    new_adj[c1, c2] = 1

        self.monomer.adjacency = new_adj

    def set_nitrogen(self, position=2):
        """
        Change an oxygen to a nitrogen.
        Example: Gal -> GalN

        Args:
            position (int): position of a carbon atom at which the bound o atom should be replaced by N

        Returns:
            rdkit id of the atom that is now a nitrogen
        """
        if position == "O":
            return

        pos = self.monomer.find_oxygen(position)
        self.monomer.structure.GetAtomWithIdx(pos).SetAtomicNum(7)
        self.monomer.x[pos, 0] = 7
        return pos

    def make_acid(self):
        """
        Make the monomer acidic by adding a acid group to the last carbon atom

        Returns:
            Nothing
        """
        emol = EditableMol(self.monomer.structure)

        c_id = int(np.argwhere(self.monomer.x[:, 1] == 6))
        o_id = EditableMol.AddAtom(emol, Atom(8))
        EditableMol.AddBond(emol, o_id, c_id, order=BondType.DOUBLE)

        self.monomer.structure = emol.GetMol()

        new_x, new_adj = self._extend_matrices(1)
        new_x[o_id, 0] = 8
        self.monomer.x = new_x

        new_adj[c_id, o_id] = 1
        new_adj[o_id, c_id] = 1
        self.monomer.adjacency = new_adj

    def desoxygenate(self, position):
        """
        desoxygenate the given position by removing the oxygen atom

        Args:
            position (int): position to desoxygenate

        Returns:
            Nothing
        """
        oxygen = self.monomer.find_oxygen(position)
        chirality = self.monomer.structure.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == position)[0])).GetChiralTag()
        print(self.monomer.structure.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == position)[0])).GetChiralTag())
        print("P1-N:", [atom.GetIdx() for atom in self.monomer.structure.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == position)[0])).GetNeighbors()])
        # self.monomer.structure.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == position)[0])).SetChiralTag(ChiralType.CHI_UNSPECIFIED)
        self._reduce_matrices([oxygen])
        print(self.monomer.structure.GetAtomWithIdx(int(np.where(self.monomer.x[:, 1] == position)[0])).GetChiralTag())
        return chirality

    def to_enantiomer(self, form):
        if self.monomer.get_isomer() == form:
            return
        ring_info = self.monomer.structure.GetRingInfo().AtomRings()[0]
        for idx in ring_info:
            if self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag() == ChiralType.CHI_TETRAHEDRAL_CCW:
                self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CW)
            elif self.monomer.structure.GetAtomWithIdx(idx).GetChiralTag() == ChiralType.CHI_TETRAHEDRAL_CW:
                self.monomer.structure.GetAtomWithIdx(idx).SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CCW)

    def _extend_matrices(self, count):
        """
        Extend the describing matrices of this monomer by count-many positions to be able to add a functional group

        Args:
            count (int): Number of atoms added in the functional group

        Returns:
            new adjacency matrix and new feature matrix
        """
        tmp_x = np.zeros((self.monomer.x.shape[0] + count, self.monomer.x.shape[1]))
        tmp_x[:self.monomer.x.shape[0], :] = self.monomer.x
        tmp_adj = np.zeros((self.monomer.adjacency.shape[0] + count, self.monomer.adjacency.shape[1] + count))
        tmp_adj[:self.monomer.adjacency.shape[0], :self.monomer.adjacency.shape[1]] = self.monomer.adjacency
        return tmp_x, tmp_adj

    def _reduce_matrices(self, positions):
        """
        Drop some atoms at the specified positions.

        Args:
            positions (List[int]): list of atom positions in the RDKit molecule to be deleted

        Returns:
            Nothing
        """
        emol = EditableMol(self.monomer.structure)

        for p in sorted(positions, reverse=True):
            self.monomer.x = np.delete(self.monomer.x, p, 0)
            self.monomer.adjacency = np.delete(self.monomer.adjacency, p, 0)
            self.monomer.adjacency = np.delete(self.monomer.adjacency, p, 1)
            EditableMol.RemoveAtom(emol, p)

        self.monomer.structure = emol.GetMol()
