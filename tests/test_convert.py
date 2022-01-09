import os.path
import sys

from glyles.converter import convert, convert_generator
from tests.test_smiles import compare_smiles


def catch_output(method, glycan=None, output_file=None, silent=True):
    logging_out, logging_err = [], []

    class writer_out(object):
        @staticmethod
        def write(data):
            logging_out.append(data)

        def close(self):
            pass

    class writer_err(object):
        @staticmethod
        def write(data):
            logging_err.append(data)

        def close(self):
            pass

    old_out = sys.stdout
    sys.stdout = writer_out()

    old_err = sys.stderr
    sys.stderr = writer_err()

    if method == convert:
        output = method(glycan=glycan, output_file=output_file, silent=silent)
    else:
        output = list(method(glycan=glycan, silent=silent))

    sys.stdout = old_out
    sys.stderr = old_err
    return output, logging_out, logging_err


def setup_test():
    def generator(x_list):
        for x in x_list:
            yield x

    def write_file(x_list, name="./test.txt"):
        f = open(name, "w")
        for x in x_list:
            f.write(x + "\n")
        f.close()
        return name

    return {
        "glycan": TestConverter.smiles_samples[0][0],
        "glycan_list": [TestConverter.smiles_samples[1][0], TestConverter.smiles_samples[2][0]],
        "glycan_file": write_file([TestConverter.smiles_samples[i][0] for i in range(3, 8)]),
        "glycan_generator": generator([TestConverter.smiles_samples[i][0] for i in range(8, 13)])
    }


def check_results(output):
    solution = [TestConverter.smiles_samples[i][1] for i in range(13)]
    i = 0
    for i, ((o_glycan, o_smiles), s_smiles) in enumerate(zip(output, solution)):
        assert o_glycan == TestConverter.smiles_samples[i][0]
        compare_smiles(o_smiles, s_smiles)
    assert i == 12

    if os.path.exists("./test.txt"):
        os.remove("./test.txt")
    if os.path.exists("./output.txt"):
        os.remove("./output.txt")


class TestConverter:
    smiles_samples = [
        ("Gal(a1-4)Gal",
         "OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),
        ("Man(a1-3)[Man(a1-6)]Man",
         "OC[C@H]3O[C@H](OC[C@H]2OC(O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H]2O)[C@@H](O)"
         "[C@@H](O)[C@@H]3O"),
        ("Man(a1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man",
         "OC[C@H]5O[C@H](O[C@H]4[C@H](O)[C@@H](CO)O[C@H](OC[C@H]3OC(O)[C@@H](O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)"
         "[C@H](O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H]3O)[C@H]4O)[C@@H](O)[C@@H](O)[C@@H]5O"),
        ("Gal(b1-4)Glc",
         "OC[C@H]2O[C@@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O"),
        ("Gal(b1-3)Glc",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@H]2O"),
        ("Glc(b1-3)Glc",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O"),
        ("Gal(b1-4)Gal",
         "OC[C@H]2O[C@@H](O[C@@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)[C@H](O)[C@@H](O)[C@H]2O"),
        ("Man(a1-4)Man",
         "OC[C@H]2O[C@H](O[C@H]1[C@H](O)[C@H](O)C(O)O[C@@H]1CO)[C@@H](O)[C@@H](O)[C@@H]2O"),
        ("Man(a1-3)Man",
         "OC[C@H]2O[C@H](O[C@@H]1[C@H](O)C(O)O[C@H](CO)[C@H]1O)[C@@H](O)[C@@H](O)[C@@H]2O"),
        ("Gal(a1-4)Gal(b1-4)Glc",
         "OC[C@H]3O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@@H](O)C(O)O[C@@H]1CO)O[C@@H]2CO)[C@H](O)"
         "[C@@H](O)[C@H]3O"),
        ("Gal(a1-3)Gal",
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@@H]1O)[C@H](O)[C@@H](O)[C@H]2O"),
        ("Glc(a1-3)Glc",
         "OC[C@H]2O[C@H](O[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O)[C@H](O)[C@@H](O)[C@@H]2O"),
        ("Man(a1-2)Man",
         "OC[C@H]2OC(O)[C@@H](O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)[C@@H]2O"),
    ]

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

    def test_convert_1_1_1(self):
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
        output, log_out, log_err = catch_output(method=convert, glycan="Glc", output_file="./invalid/path/out.txt", silent=True)

        assert output is None
        assert len(log_out) == 4
        assert log_out[0] == "Glc"
        compare_smiles(log_out[2], "[C@H]1(O)[C@H](O)[C@@H](O)C(O)O[C@@H]1CO")
        assert len(log_err) == 0

    def test_convert_2_2(self):
        output, log_out, log_err = catch_output(method=convert, glycan="Glc", output_file="./invalid/path/out.txt", silent=False)

        assert output is None
        assert len(log_out) == 4
        assert log_out[0] == "Glc"
        compare_smiles(log_out[2], "[C@H]1(O)[C@H](O)[C@@H](O)C(O)O[C@@H]1CO")
        assert len(log_err) == 2
        assert log_err[0] == "Path of output-file does not exist! Results will be printed on stdout."

    def test_convert_3_1(self):
        output, log_out, log_err = catch_output(method=convert, glycan="Man", silent=True)

        assert output is None
        assert len(log_out) == 4
        assert log_out[0] == "Man"
        assert len(log_err) == 0

    def test_convert_3_2(self):
        output, log_out, log_err = catch_output(method=convert, glycan="Man", silent=False)

        assert output is None
        assert len(log_out) == 6
        assert log_out[0] == "No output-file specified, results will be printed on stdout."
        assert log_out[2] == "Man"
        assert len(log_err) == 0

