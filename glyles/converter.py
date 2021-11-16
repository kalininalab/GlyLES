import os
import sys

from glyles.grammar.parse import Glycan


def convert(glycan=None, glycan_list=None, glycan_file=None, output_file=None, generator=False, silent=True):
    """
    General user interaction interface to use this library.

    idea: Implement a generator argument that arbitrarily generates a stream of glycans

    Args:
        glycan (str): Single glycan to be converted from IUPAC to SMILES
        glycan_list (List[str]): list of glycans to convert
        glycan_file (str): File to read the glycans from
        output_file (str): File to save the conversions in
        generator (bool): Flag indicating the use of this method as generator
        silent (bool): Flag indicating to have no output from this method

    Returns:
        Nothing
    """
    glycans = []

    # fill the list of glycans to convert
    if glycan is not None:
        glycans.append(glycan)
    if glycan_list is not None:
        glycans += glycan_list
    if glycan_file is not None:
        if not os.path.isfile(glycan_file):
            pass
        for line in open(glycan_file, "r").readlines():
            glycans.append(line.strip())
    if len(glycans) == 0:
        if not silent:
            print("List of glycans is empty")
        return

    # determine the output format
    # TODO: Check for generator, a generator makes this unnecessary
    if output_file is not None:
        if os.path.isdir(os.path.dirname(output_file)):
            output = open(output_file, "w")
        else:
            if not silent:
                print("Path of output-file does not exist! Results will be printed on stdout.", file=sys.stderr)
            output = sys.stdout
    else:
        if not silent:
            print("No output-file specified, results will be printed on stdout.")
        output = sys.stdout

    # Convert the glycans ...
    for glycan in glycans:
        try:
            # ... by passing them to the glycan class to parse them ...
            smiles = Glycan(glycan).get_smiles()

            # ... and return them as intended
            if generator:
                yield smiles
            else:
                print(glycan, smiles, file=output, sep=",")
        # catch any exception at glycan level to not destroy the whole pipeline because of one mis-formed glycan
        except Exception as e:
            print(f"An exception occurred with {glycan}:", e.__class__, file=sys.stderr)
            print("Error message:", e.__str__(), file=sys.stderr)
            if generator:
                yield ""

    output.close()


if __name__ == '__main__':  # testing
    for smiles in convert(glycan="Man(a1-2)Man", output_file="./blbabl/out.txt"):
        print(smiles)
