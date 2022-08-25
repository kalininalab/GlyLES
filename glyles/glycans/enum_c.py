import numpy as np

from glyles.glycans.utils import UnreachableError, Tree


def enumerate_carbon(monomer):
    """
    Enumerate carbon atoms starting with the one that are also part of the root monomer and then their side chains in
    order of the original carbon atoms.

    Args:
        monomer (Monomer): object of class monomer to enumerate the carbon atoms for

    Returns:
        Nothing
    """
    # first enumerate all ring carbons and the attachments at the first and last ring atom
    c_atoms = np.where((monomer.x[:, 0] == 6) & (monomer.x[:, 2] == 1) & (monomer.x[:, 3] == 1))[0]
    ring_o = np.where((monomer.x[:, 0] == 8) & (monomer.x[:, 2] == 1) & (monomer.x[:, 3] == 1))[0]
    if ring_o.size == 0:
        next_c_id = enumerate_c_atoms(monomer, c_atoms, -1)
    else:
        next_c_id = enumerate_c_atoms(monomer, c_atoms, ring_o)

    # then iterate over all those carbons and enumerate them beginning at C1
    for c_id in range(1, next_c_id):
        c_index = int(np.where(monomer.x[:, 1] == c_id)[0])
        candidates = np.where((monomer.adjacency[c_index, :] != 0) & (monomer.x[:, 1] == 0))[0]
        if candidates.size != 0:
            for candidate in list(candidates):
                if sum(monomer.adjacency[candidate, :]) > 1:
                    next_c_id = enumerate_side_chain(monomer, c_index, candidate, next_c_id)


def enumerate_side_chain(monomer, parent, atom, next_c_id):
    """
    Not enumerating ring atoms

    Args:
        monomer (Monomer): object of class monomer to enumerate the carbon atoms for
        parent (int): RDKit ID of the parent atom of the one to look at here
        atom (int): RDKit ID of the atom to look at in this recursive call
        next_c_id (int): ID to assign to the next carbon atom

    Returns:
        Number (not RDKit ID) of the highest carbon atom in the molecule
    """
    # if side chain starts with a carbon, put a number on it
    if monomer.x[atom, 0] == 6:
        monomer.x[atom, 1] = next_c_id
        next_c_id += 1

    # identify children
    candidates = np.where((monomer.adjacency[atom, :] != 0) & (monomer.x[:, 1] == 0))[0].tolist()
    if parent in candidates:
        candidates.remove(parent)

    # and iterate over them to enumerate them in order as they appear
    if len(candidates) > 0:
        for candidate in candidates:
            next_c_id = enumerate_side_chain(monomer, atom, candidate, next_c_id)

    # return the ID of the next carbon that we see
    return next_c_id


def enumerate_c_atoms(monomer, c_atoms, ringo):
    """
    Enumerate all carbon atoms starting from the first one

    Args:
        monomer (Monomer): object of class monomer to enumerate the carbon atoms for
        c_atoms (List[int]): List of all ids of C atoms in the ring
        ringo (int): id of the oxygen atom in the ring of the monomer

    Returns:
        Number (not RDKit ID) of the highest carbon atom in the molecule
    """
    if len(c_atoms) > 0:
        # create a tree of all carbon atoms directly connected to the main ring of the monomer
        c_tree = Tree()
        stack = [(-1, c_atoms[0])]
        while len(stack) != 0:
            p_id, c_id = stack[-1]
            stack = stack[:-1]
            c_tree.add_node(c_id, p_id)

            children = np.where((monomer.adjacency[c_id, :] == 1) & (monomer.x[:, 0] == 6) & (monomer.x[:, 3] == 1))[0]
            for c in children:
                if int(c) not in c_tree.nodes:
                    stack.append((c_id, int(c)))

        # find the deepest node and rehang the tree to this node
        deepest_id, _ = c_tree.deepest_node()
        c_tree = c_tree.rehang_tree(deepest_id)
        longest_c_chain = c_tree.longest_chain()
        assert len(longest_c_chain) >= len(c_atoms), "Something went wrong, AGAIN"

        # now the two C1 candidates can be found at the ends of the longest chain
        start, end = longest_c_chain[0], longest_c_chain[-1]

        # check conditions
        start_o_conn, end_o_conn = \
            np.argwhere((monomer.adjacency[start, :] == 1) & (monomer.x[:, 0] == 8) & (monomer.x[:, 2] != 1) &
                        (monomer.x[:, 3] == 1)).squeeze().size > 0, \
            np.argwhere((monomer.adjacency[end, :] == 1) & (monomer.x[:, 0] == 8) & (monomer.x[:, 2] != 1) &
                        (monomer.x[:, 3] == 1)).squeeze().size > 0

        # decide on c1
        if start_o_conn and end_o_conn:
            if not evaluate_distance(monomer, start, end, ringo):
                longest_c_chain = reversed(longest_c_chain)
        elif end_o_conn:
            longest_c_chain = reversed(longest_c_chain)
    else:
        longest_c_chain = monomer.c1_find(monomer.structure)
    longest_c_chain = list(longest_c_chain)

    # enumerate along chain
    next_c_id = 1
    for c in longest_c_chain:
        monomer.x[c, 1] = next_c_id
        next_c_id += 1

    return next_c_id


def evaluate_distance(monomer, start, end, ringo):
    """
    Try to decide on C1 based on their distance to the oxygen in the ring

    Args:
        monomer (Monomer): object of class monomer to enumerate the carbon atoms for
        start (int): id of the first candidate for C1 atom
        end (int): id of the second candidate for C1 atom
        ringo (int): index of the oxygen atom in the ring

    Returns:
        Bool indicating that the start id is the C1 atom
    """
    adj = monomer.adjacency.copy()

    # as we have an adjacency matrix, multiply it with itself until one of the fields is non-zero
    while adj[start, ringo] == 0 and adj[end, ringo] == 0:
        adj = adj @ monomer.adjacency

    # if both fields are non-zero, we cannot decide here and have to go further
    if adj[start, ringo] > 0 and adj[end, ringo] > 0:
        return equidistant(monomer, start, end)

    return adj[start, ringo] > 0


def equidistant(monomer, start, end):
    """
    Decider for C1 in case the previous splitting rules were all tied.

    Args:
        monomer (Monomer): object of class monomer to enumerate the carbon atoms for
        start (int): id of the first candidate for C1 atom
        end (int): id of the second candidate for C1 atom

    Returns:
        Bool indicating that the start id is the C1 atom
    """
    # determine first carbon atom in the ring
    c_start_candidates = np.where((monomer.adjacency[start, :] == 1) &
                                  (monomer.x[:, 0] == 6) & (monomer.x[:, 2] == 1))[0]
    c_end_candidates = np.where((monomer.adjacency[end, :] == 1) & (monomer.x[:, 0] == 6) & (monomer.x[:, 2] == 1))[0]

    if c_start_candidates.size == 1 and c_end_candidates.size == 1:
        start_ring_c = int(c_start_candidates)
        end_ring_c = int(c_end_candidates)

        # check if those ring carbons have an attached oxygen
        start_ring_c_o_candidates = np.where((monomer.adjacency[start_ring_c, :] == 1) &
                                             (monomer.x[:, 0] == 8) & (monomer.x[:, 2] != 1))[0]
        end_ring_c_o_candidates = np.where((monomer.adjacency[end_ring_c, :] == 1) &
                                           (monomer.x[:, 0] == 8) & (monomer.x[:, 2] != 1))[0]

        if start_ring_c_o_candidates.size == 1 and end_ring_c_o_candidates.size == 1:
            raise UnreachableError("C1 atom cannot be detected")
        elif start_ring_c_o_candidates.size == 1:
            return True
    return c_start_candidates.size == 1
