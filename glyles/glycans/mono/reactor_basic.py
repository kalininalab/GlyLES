import copy
import re

import numpy as np

from glyles.glycans.factory.factory_o import OpenFactory
from glyles.grammar.GlycanLexer import GlycanLexer


def check_for_resizing(monomer, names, types):
    """
    Change the monosaccharide in case it's indicated that the monosaccharide has more carbon atoms than normal.

    Args:
        monomer (Monomer): Monomer to be modified according to functional groups
        names (List[str]): List of the names of the functional groups
        types (List[int]): List of type definitions from the recipe of a monosaccharide

    Returns:
        Nothing
    """
    sac_index = types.index(GlycanLexer.SAC)

    # if there's nothing to do, return
    if len(types) == sac_index + 1 or (len(types) > sac_index + 1 and types[sac_index + 1] != GlycanLexer.SAC):
        return

    # identify the root and the new size of the monomer by identifying the indices of their descriptions
    if names[sac_index + 1] in ["Pen", "Hex", "Hep", "Oct"]:
        len_index = sac_index + 1
    else:
        len_index = sac_index
        sac_index = len_index + 1

    # extract the orientations of the hydroxy-groups along the enlarged chain
    if sac_index != 0 and names[sac_index - 1].count("L") + \
            names[sac_index - 1].count("D") == len(names[sac_index - 1]):
        orientations = names[sac_index - 1]
    else:
        orientations = ""

    # determine the chain that has to be added to the last carbon atom
    c_count = np.count_nonzero(monomer.x[:, 0] == 6)
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

    # insert the orientations
    for c in orientations[:-1]:
        if c == "L":
            extension = extension.replace("[C?H]", "[C@@H]", 1)
        if c == "D":
            extension = extension.replace("[C?H]", "[C@H]", 1)
    extension = extension.replace("[C?H]", "C")

    # find oxygen of c6 to delete it and find c6 to replace it
    c_id = int(np.where(monomer.x[:, 1] == int(max(monomer.x[monomer.x[:, 0] == 6, 1])))[0])
    # if the carbon to be extended is part of the ring, put the extension into a side_chain
    if monomer.x[c_id, 2] & 0b1:
        extension = extension.replace(")", ")(", 1) + ")"
    # if the original molecule has no oxygen at its last carbon, remove the first oxygen from the extension
    if monomer.x[((monomer.adjacency[c_id] - monomer.x[:, 2]) > 0) & (monomer.x[:, 0] == 8), 0].size == 0:
        extension = extension.replace("(O)", "", 1)
    monomer.structure.GetAtomWithIdx(monomer.find_oxygen(c_count)).SetAtomicNum(50)
    monomer.structure.GetAtomWithIdx(c_id).SetAtomicNum(32)

    # generate SMILES from rings oxygen and replace c6's oxygen and c6 by the extension
    smiles = monomer.to_smiles(0, root_idx=10)
    smiles = smiles.replace("[SnH]", "").replace("[GeH2]", extension).replace("[GeH]", extension).replace("()", "")

    # update the monomer
    monomer.smiles = smiles
    monomer.structure = None
    monomer.get_structure()


def check_for_open_form(monomer, names, types, longer=False):
    """
    This method checks if the monomer is converted into an open form.

    Args:
        monomer (Monomer): Monomer to be modified according to functional groups
        names (List[str]): List of the names of the functional groups
        types (List[int]): List of type definitions from the recipe of a monosaccharide
        longer (bool): True, if the glycan is longer/shorted than normal

    Returns:
        Nothing
    """
    # identify the descriptions of what to do and determine the new SMILES string
    sac_index = types.index(GlycanLexer.SAC)
    params = copy.copy(OpenFactory()[names[sac_index] + "-ol"])

    if longer:
        last = params["smiles"].rfind("C")
        count = len(re.findall(r"\(C\)", params["smiles"]))
        if "Pen" in names:
            maximum = 5
        elif "Hex" in names:
            maximum = 6
        elif "Hep" in names:
            maximum = 7
        else:  # if "Oct" in names
            maximum = 8
        params["smiles"] = params["smiles"][:last] + "C(O)" * max((maximum - count), 0) + params["smiles"][last:]

    if "-onic" in names or "-aric" in names or "-ulosaric":
        params["smiles"] = params["smiles"].replace("C", "C(=O)", 1)
    if "-aric" in names or "-ulosaric" in names and names[sac_index] != "Qui":
        params["smiles"] = params["smiles"][:-1] + "(=O)O"
    if "-ulosonic" in names or "-olusaric" in names:
        params["smiles"] = re.sub(r"OC\[C@*H]\(O\)", "OC(=O)C(=O)", params["smiles"], count=1)

    # update the monomer accordingly
    monomer.smiles = params["smiles"]
    monomer.c1_find = params["c1_find"]
    monomer.structure = None
    monomer.get_structure()


def change_base_monomer(monomer, names, types):
    """
    Check if the basic structure of the monomer has been changed

    Args:
        monomer (Monomer): Monomer to be modified according to functional groups
        names (List[str]): List of the names of the functional groups
        types (List[int]): List of type definitions from the recipe of a monosaccharide

    Returns:
        modified monomer object
    """
    longer = len(set(names).intersection({"Pen", "Hex", "Hep", "Oct"})) != 0
    open_form = len(set(names).intersection({"-ol", "-onic", "-aric", "-ulosonic", "-ulosaric"})) != 0

    if not longer and not open_form:
        # nothing to do
        return monomer

    if longer and not open_form:
        # parse for sth. like LDManHep or DDManHep -> afterwards, no need to look for LD/DD/... and Hep/Hex/Pen/Oct
        check_for_resizing(monomer, names, types)
        return monomer

    # open form glycan, changed ends and eventually changed length
    check_for_open_form(monomer, names, types, longer)

    return monomer
