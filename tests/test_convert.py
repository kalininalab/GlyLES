import os

from glyles.converter import convert, convert_generator
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


class TestConverter:
    def test_file_output(self):
        args = setup_test()

        convert(args["glycan"], args["glycan_list"], args["glycan_file"], args["glycan_generator"],
                output_file="./output.txt")
        output = [line.split(",") for line in open("./output.txt", "r").readlines()]
        check_results(output)

    def test_return_output(self):
        args = setup_test()

        output = convert(args["glycan"], args["glycan_list"], args["glycan_file"], args["glycan_generator"],
                         returning=True)
        check_results(output)

    def test_generator(self):
        args = setup_test()

        output = convert_generator(args["glycan"], args["glycan_list"], args["glycan_file"], args["glycan_generator"])
        check_results(output)

    def test_orientation_1(self):
        output = convert("Glc a", returning=True)

        assert len(output) == 1
        assert output[0][0] == "Glc a"
        compare_smiles(output[0][1], "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O")

    def test_orientation_2(self):
        output = convert("Glc b", returning=True)

        assert len(output) == 1
        assert output[0][0] == "Glc b"
        compare_smiles(output[0][1], "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O")

    """def test_convert_1_1_1(self):
        output, log_out, log_err = catch_output(method=convert, silent=False)

        assert output is None
        assert len(log_out) == 0
        assert len(log_err) == 2

    def test_convert_1_1_2(self):
        output, log_out, log_err = catch_output(method=convert, silent=False)

        assert output is None
        assert len(log_out) == 0
        assert len(log_err) == 2
        assert log_err[0] == "List of glycans is empty"

    def test_convert_1_2_1(self):
        output, log_out, log_err = catch_output(method=convert_generator, silent=True)

        assert len(output) == 0
        assert len(log_out) == 0
        assert len(log_err) == 0

    def test_convert_1_2_2(self):
        output, log_out, log_err = catch_output(method=convert_generator, silent=False)

        assert len(output) == 0
        assert len(log_out) == 0
        assert len(log_err) == 2
        assert log_err[0] == "List of glycans is empty"

    def test_convert_2_1(self):
        output, log_out, log_err = catch_output(method=convert, glycan="Glc", output_file="./invalid/path/out.txt",
                                                silent=True)

        assert output is None
        assert len(log_out) == 4
        assert log_out[0] == "Glc"
        compare_smiles(log_out[2], "[C@H]1(O)[C@H](O)[C@@H](O)C(O)O[C@@H]1CO")
        assert len(log_err) == 0

    def test_convert_2_2(self):
        output, log_out, log_err = catch_output(method=convert, glycan="Glc", output_file="./invalid/path/out.txt",
                                                silent=False)

        assert output is None
        assert len(log_out) == 4
        assert log_out[0] == "Glc"
        compare_smiles(log_out[2], "[C@H]1(O)[C@H](O)[C@@H](O)C(O)O[C@@H]1CO")
        assert len(log_err) == 2
        assert log_err[0] == "Path of output-file does not exist! Results will be printed on stdout."

    def test_convert_3_1(self):
        output, log_out, log_err = catch_output(method=convert, glycan="Man", silent=True, returning=False)

        assert output is None
        assert len(log_out) == 4
        assert log_out[0] == "Man"
        assert len(log_err) == 0

    def test_convert_3_2(self):
        output, log_out, log_err = catch_output(method=convert, glycan="Man", silent=False, returning=False)

        assert output is None
        assert len(log_out) == 6
        assert log_out[0] == "No output-file specified, results will be printed on stdout."
        assert log_out[2] == "Man"
        assert len(log_err) == 0"""
