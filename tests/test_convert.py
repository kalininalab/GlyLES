import os
import time

import numpy as np
import pytest
from PIL import Image

from glyles import Glycan
from glyles.converter import convert, convert_generator
from glyles.glycans.utils import sanitize_smiles
from tests.utils import setup_test, smiles_samples
from rdkit import Chem


def check_results(output):
    solution = [smiles_samples[i][1] for i in range(13)]
    i = 0
    for i, ((o_glycan, o_smiles), s_smiles) in enumerate(zip(output, solution)):
        assert o_glycan == smiles_samples[i][0]
        compare_smiles(o_smiles, s_smiles)
    assert i == 12

    if os.path.exists("./test.txt"):
        os.remove("./test.txt")
    if os.path.exists("./output.txt"):
        os.remove("./output.txt")


def compare_smiles(computed, solution):
    c = Chem.MolFromSmiles(computed)
    Chem.Kekulize(c)
    c_rdkit = Chem.MolToSmiles(c, kekuleSmiles=True)

    s = Chem.MolFromSmiles(solution)
    Chem.Kekulize(s)
    s_rdkit = Chem.MolToSmiles(s, kekuleSmiles=True)

    assert c_rdkit == s_rdkit


@pytest.mark.parametrize("iupac", [
    "Man",
    "Man(a1-4)Man",
    "Man(a1-4)[Man(a1-3)]Man",
    "Man(a1-4)[Man(a1-3)][Man(a1-2)]Man",
])
def test_demo(iupac):
    output = convert(iupac)

    assert len(output) == 1
    assert len(output[0]) == 2
    assert output[0][0] == iupac
    assert all(c in output[0][1] for c in {"C", "O", "[", "]", "@"})

@pytest.mark.parametrize("iupac", [
    "Man(?1-4)Man",
    "Man(a1-?)Man",
    "Man(?1-?)Man",
    "{Man(a1-4)}Man(a1-2)Man",
    "{Man(?1-4)}Man(a1-?)Man",
    "{Man(?1-?)}Man(?1-?)Man",
])
def test_demo_no_resolution(iupac):
    output = convert(iupac)

    assert len(output) == 1
    assert len(output[0]) == 2
    assert output[0][0] == iupac
    assert output[0][1] == ""

def test_file_output():
    args = setup_test()

    convert(args["glycan"], args["glycan_list"], args["glycan_file"], args["glycan_generator"],
            output_file="./output.txt")
    output = [line.split(",") for line in open("./output.txt", "r").readlines()]
    check_results(output)

def test_return_output():
    args = setup_test()

    output = convert(args["glycan"], args["glycan_list"], args["glycan_file"], args["glycan_generator"],
                        returning=True)
    check_results(output)

def test_generator():
    args = setup_test()

    output = convert_generator(args["glycan"], args["glycan_list"], args["glycan_file"], args["glycan_generator"])
    check_results(output)

def test_orientation_1():
    output = convert("Glc a", returning=True)

    assert len(output) == 1
    assert output[0][0] == "Glc a"
    compare_smiles(output[0][1], "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O")

def test_orientation_2():
    output = convert("Glc b", returning=True)

    assert len(output) == 1
    assert output[0][0] == "Glc b"
    compare_smiles(output[0][1], "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O")

@pytest.mark.slow
def test_parallel_list_input():
    with open("data/runtime.tsv", "r") as runtime:
        glycans = [line.strip().split("\t")[0] for line in runtime.readlines()]
    start = time.time()
    smiles1 = convert(glycan_list=glycans)
    mid = time.time()
    smiles2 = convert(glycan_list=glycans, cpu_count=14)
    end = time.time()

    assert len(smiles1) == len(smiles2) == len(glycans)
    assert (end - mid) < (mid - start)

@pytest.mark.slow
def test_parallel_gen_input():
    with open("data/runtime.tsv", "r") as runtime:
        glycans = [line.strip().split("\t")[0] for line in runtime.readlines()]
    start = time.time()
    smiles1 = convert(glycan_generator=iter(glycans))
    mid = time.time()
    smiles2 = convert(glycan_generator=iter(glycans), cpu_count=14)
    end = time.time()

    assert len(smiles1) == len(smiles2) == len(glycans)
    assert (end - mid) < (mid - start)

def test_smiles_clean():
    assert sanitize_smiles("SDJCBPIOUCODJCOBC") == "SDJCBPIOUCODJCOBC"
    assert sanitize_smiles("DIC(DONC)WOUC") == "DIC(DONC)WOUC"
    assert sanitize_smiles("DPIUCDBPSIDU((CPIDBC)PID)") == "DPIUCDBPSIDU(CPIDBCPID)"
    assert sanitize_smiles("SDJC((PSODUCBN))SOD:C") == "SDJC(PSODUCBN)SOD:C"
    assert sanitize_smiles("A:DO(C(OPDIC))PODUC") == "A:DO(COPDIC)PODUC"

def test_summary():
    output = Glycan("Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-3)[Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-6)]"
                    "Gal(b1-3)[GlcNAc(a1-4)Gal(b1-4)GlcNAc6S(b1-6)]GalNAc", tree_only=True).summary()
    assert output["formula"] == "C92H153N7O67S"
    assert 2459 < output["weight"] < 2460
    assert output["atoms"] == 167
    assert output["bonds"] == 179
    assert output["rings"] == 13
    assert output["monomers"] == 13
    assert output["types"]["GalNAc"] == 3
    assert output["types"]["GlcNAc6S"] == 1
    assert output["types"]["Gal"] == 4
    assert output["types"]["GlcNAc"] == 3
    assert output["types"]["Fuc"] == 2
    assert output["root"] == "GalNAc"
    assert output["leaves"].count("GlcNAc") == 1
    assert output["leaves"].count("GalNAc") == 2
    assert output["leaves"].count("Fuc") == 2
    assert output["depth"] == 4
