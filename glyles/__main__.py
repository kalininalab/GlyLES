import os
import sys
import argparse

from glyles import convert
from glyles.version import __version__


def parse_args(args):
    """
    Parse the arguments given to this script.

    Args:
        args: Arguments taken from sys.argv, provided explicitly to be able to test the functionality of this method

    Returns:
        The arguments parsed from the commandline into a dictionary as in kwargs
    """
    parser = argparse.ArgumentParser(
        prog="glyles",
        description="A tool to convert IUPAC representation of Glycans into SMILES representation"
    )
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        action="store",
        help="Input IUPAC string, list or file path",
        nargs='+'
    )
    parser.add_argument(
        '-o',
        '--output',
        required=True,
        action="store",
        help="Output file path"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )
    return vars(parser.parse_args(args))


def parse_list(input_list):
    """
    Parse a list of different inputs into a list.

    Args:
        input_list (List[str]): List of IUPAC-condensed strings and filenames to files containing the glycan names

    Returns:
        A list with glycan names taken from the files and the individual glycan names
    """
    glycans = []
    for x in input_list:
        if os.path.isfile(x):
            with open(x, "r") as data:
                glycans += [x.strip() for x in data.readlines()]
        else:
            glycans.append(x)
    return glycans


def main(args):
    """
    Main functionality for the CLI of GlyLES.

    Args:
        args: Arguments taken from sys.argv, provided explicitly to be able to test the functionality of this method
    """
    # parse the arguments
    kwargs = parse_args(args)

    # check if the output already exists
    if os.path.isfile(kwargs["output"]):
        if input(
                "The given output file already exist, do you want to proceed and overwrite it? [y]/n\n"
        ) not in ["", "Y", "y"]:
            exit(0)

    # check if it's a one-element list, in case replace the input with the only element in the input list
    if isinstance(kwargs["input"], list) and len(kwargs["input"]) == 1:
        kwargs["input"] = kwargs["input"][0]

    # forward the input to glyles ...
    if isinstance(kwargs["input"], list):
        data = parse_list(kwargs["input"])
        convert(glycan_list=data, output_file=kwargs["output"])
    elif os.path.isfile(kwargs["input"]):
        convert(glycan_file=kwargs["input"], output_file=kwargs["output"])
    elif isinstance(kwargs["input"], str):
        convert(glycan=kwargs["input"], output_file=kwargs["output"])
    else:
        # ... or raise an error if no matching input type has been found
        print(
            f"An error occurred with the input type: {type(kwargs['input'])}.\n"
            f"Please open an issue on github for this.\n"
            f"https://github.com/kalininalab/GlyLES/issues",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main(sys.argv[1:])
