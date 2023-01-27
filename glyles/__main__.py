import os
import sys
import argparse

# SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# sys.path.append(os.path.dirname(SCRIPT_DIR))

from glyles import convert


def main():
    input_string_flag = False
    input_list_flag = False
    input_path_flag = False

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
        '--override',
        required=False,
        action="store_true",
        help="Override output file"
    )
    # parser.add_argument("-h", "--help", action="help")
    # parser.add_argument("-v", "--version", action="version")
    args = parser.parse_args()

    input = args.input
    output = args.output

    if len(input) == 1:
        input = input[0]

    if os.path.isfile(input):
        input_path_flag = True
    else:
        if type(input) == list:
            input_list_flag = True

        elif type(input) == str:
            input_string_flag = True

    if os.path.isfile(output):
        if not args.override:
            print("Given output exist, please pick another name or use --override")
            sys.exit()

    if input_string_flag:
        convert(glycan=input, output_file=output)

    elif input_list_flag:
        convert(glycan_list=input, output_file=output)

    elif input_path_flag:
        convert(glycan_file=input, output_file=output)


if __name__ == "__main__":
    main()
