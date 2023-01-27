import os
import sys
import argparse
import pkg_resources


from glyles import convert


def main():
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
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {pkg_resources.get_distribution('glyles').version}"
    )
    args = parser.parse_args()
    if os.path.isfile(args.output):
        if not args.override:
            print("Given output exist, please pick another name or use --override")
            sys.exit()

    if len(args.input) == 1:
        args.input = args.input[0]

    if os.path.isfile(args.input):
        convert(glycan_file=args.input, output_file=args.output)
    elif isinstance(args.input, list):
        convert(glycan_list=args.input, output_file=args.output)
    elif type(input) == str:
        convert(glycan=args.input, output_file=args.output)


if __name__ == "__main__":
    main()
