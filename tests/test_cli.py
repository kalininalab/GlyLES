import os

from rdkit import Chem

from glyles.__main__ import main
from tests.test_derivatives import valid_atomic_nums


def compare_smiles(computed, solution):
    c = Chem.MolFromSmiles(computed)
    assert all([a.GetAtomicNum() in valid_atomic_nums for a in c.GetAtoms()])
    Chem.Kekulize(c)
    c_rdkit = Chem.MolToSmiles(c, kekuleSmiles=True)

    s = Chem.MolFromSmiles(solution)
    Chem.Kekulize(s)
    s_rdkit = Chem.MolToSmiles(s, kekuleSmiles=True)

    assert c_rdkit == s_rdkit


def test_cli():
    main([
        # "glyles",
        "-i",
        "data/test_cli1.txt",
        "ManA",
        "Glc4P",
        "data/test_cli2.txt",
        "data/test_cli3.txt",
        "Fru2Pfb",
        "Man6Pb",
        "-o",
        "data/test_cli_output.txt",
    ])

    assert os.path.isfile("data/test_cli_output.txt")
    with open("data/test_cli_output.txt", "r") as data:
        output = [tuple(line.strip().split(",")) for line in data.readlines()]

    with open("data/test_cli_solution.txt", "r") as data:
        solution = [tuple(line.strip().split("\t")) for line in data.readlines()]

    assert len(solution) == len(output)
    for (sol_iupac, sol_smiles), (out_iupac, out_smiles) in zip(solution, output):
        assert sol_iupac == out_iupac
        compare_smiles(out_smiles, sol_smiles)
