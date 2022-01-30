import numpy as np
from rdkit.Chem.rdchem import Atom, BondType, EditableMol

from glyles.grammar.GlycanLexer import GlycanLexer


def not_implemented_message(mod):
    print(f"ModificationNotImplementedWarning: {mod} Modification not implemented. The returned molecule will not have "
          f"this modification")


class Reactor:
    """
    Class with access to protected classes fields managing the modifications of a root monomer.
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

        Args:
            names (List[str]): name (string representation) of the modification
            types (List[int]): Type of the parsed stings based on GlycanLexer.TYPE

        Returns:
            Modified monomer
        """
        for name, t in zip(names, types):
            if t != GlycanLexer.MOD:
                continue
            if len(name) == 1:
                if name[0] == "N":
                    self.set_nitrogen()
                if name[0] == "A":
                    self.make_acid()
            elif len(name) == 2:
                if name[0].isdigit():
                    if name[1] == "d":  # ?d
                        not_implemented_message(name)
                    elif name[1] == "e":  # ?d
                        not_implemented_message(name)
                    elif name[1] == "S":
                        self.add_sulfur(int(name[0]))
                    elif name[1] == "P":
                        self.add_phosphate(int(name[0]))
                    elif name[1] == "N":
                        self.set_nitrogen(position=int(name[0]))
                else:
                    if name == "Ac":
                        self.add_acid(position=5)
            elif len(name) == 3:
                if name == "NAc":
                    self.add_acid(pos=self.set_nitrogen())
                elif name[0].isdigit():
                    if name.endswith("Me"):
                        self.add_methyl(position=int(name[0]))
                    elif name.endswith("Ac"):
                        self.add_acid(position=int(name[0]))
            elif len(name) == 7:
                if name.endswith("Me-"):
                    self.add_methyl(position=int(name[0]))
                if name.endswith("Ac-"):
                    self.add_acid(position=int(name[0]))

        return self.monomer

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
            pos = self.monomer.find_oxygen(position, check_for=[8, 7])

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

    def add_methyl(self, position):
        """

        Args:
            position:

        Returns:

        """
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

    def set_nitrogen(self, position=2):
        """
        Change an oxygen to a nitrogen.
        Example: Gal -> GalN

        Args:
            position (int): position of a carbon atom at which the bound o atom should be replaced by N

        Returns:
            rdkit id of the atom that is now a nitrogen
        """
        pos = self.monomer.find_oxygen(position, check_for=[8, 7])
        self.monomer.structure.GetAtomWithIdx(pos).SetAtomicNum(7)
        self.monomer.x[pos, 0] = 7
        return pos

    def make_acid(self):
        """

        Returns:

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

    def _extend_matrices(self, count):
        """

        Args:
            count:

        Returns:

        """
        tmp_x = np.zeros((self.monomer.x.shape[0] + count, self.monomer.x.shape[1]))
        tmp_x[:self.monomer.x.shape[0], :] = self.monomer.x
        tmp_adj = np.zeros((self.monomer.adjacency.shape[0] + count, self.monomer.adjacency.shape[1] + count))
        tmp_adj[:self.monomer.adjacency.shape[0], :self.monomer.adjacency.shape[1]] = self.monomer.adjacency
        return tmp_x, tmp_adj