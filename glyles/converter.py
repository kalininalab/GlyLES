import os
import sys
import logging

from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.utils import ParseError
from glyles.grammar.parse import Glycan


def preprocess_glycans(glycan, glycan_list, glycan_file):
    """
    Preprocess the static inputs for the parsing into one single list

    Args:
        glycan (str): single glycan to parse
        glycan_list (List[str]): list of glycans to parse
        glycan_file (str): filepath of file to read glycans from

    Returns:
        list of glycans in the order they are handed to the function, i.e. glycan, glycan_list, glycan_file
    """
    glycans = []

    # fill a list with all glycans to convert
    if glycan is not None:
        glycans.append(glycan)
    if glycan_list is not None:
        glycans += glycan_list
    if glycan_file is not None:
        # check if the file is valid and read it out
        if not os.path.isfile(glycan_file):
            raise ValueError(f"{glycan_file} does not exists, cannot read glycans.")
        for line in open(glycan_file, "r").readlines():
            glycans.append(line.strip())
    return glycans


def convert(
        glycan=None,
        glycan_list=None,
        glycan_file=None,
        glycan_generator=None,
        output_file=None,
        returning=True,
        verbose=logging.INFO,
        full=True,
):
    """
    General user interaction interface to use this library.

    Args:
        glycan (str): Single glycan to be converted from IUPAC to SMILES
        glycan_list (List[str]): list of glycans to convert
        glycan_file (str): File to read the glycans from
        glycan_generator (generator): generator yielding iupac representation.
            Together with output_generator=True this does not create any lists
        output_file (str): File to save the converted glycans in
        returning (bool): Flag indicating to return a list of tuples
        verbose (int): Flag indicating to have no prints from this method
        full (bool): Flag indicating that only fully convertible glycans should be returned, i.e. all modifications
            such as 3-Anhydro-[...] are also present in the SMILES

    Returns:
        List of type (IUPAC, SMILES) items giving the converted SMILES formulas. Only if returning=True is set.
    """

    logging.basicConfig(level=verbose)

    # collect all data and return if no data were provided
    glycans = preprocess_glycans(glycan, glycan_list, glycan_file)
    if len(glycans) == 0 and glycan_generator is None:
        logging.info("List of glycans is empty")
        return

    # determine the output format
    if output_file is not None:
        if os.path.isdir(os.path.dirname(os.path.abspath(output_file))):
            output = open(output_file, "w")
        else:
            logging.warning("Path of output-file does not exist! Results will be printed on stdout.")
            output = sys.stdout
        returning = False
    else:
        if returning:
            output = []
        else:
            logging.warning("No output-file specified, results will be printed on stdout.")
            output = sys.stdout

    # convert the IUPAC strings into SMILES strings from the input list
    if len(glycans) != 0:
        for glycan, smiles in convert_generator(glycan_list=glycans, verbose=verbose, full=full):
            if returning:
                output.append((glycan, smiles))
            else:
                print(glycan, smiles, file=output, sep=",")

    # and from the input generator
    if glycan_generator is not None:
        for glycan, smiles in convert_generator(glycan_generator=glycan_generator, verbose=verbose, full=full):
            if returning:
                output.append((glycan, smiles))
            else:
                print(glycan, smiles, file=output, sep=",")

    if returning:
        return output
    elif output_file is not None:
        output.close()


def convert_generator(
        glycan=None,
        glycan_list=None,
        glycan_file=None,
        glycan_generator=None,
        verbose=logging.INFO,
        full=True,
):
    """
    General user interaction interface to use this library.

    Args:
        glycan (str): Single glycan to be converted from IUPAC to SMILES
        glycan_list (List[str]): list of glycans to convert
        glycan_file (str): File to read the glycans from
        glycan_generator (generator): generator yielding iupac representation.
            Together with output_generator=True this does not create any lists
        verbose (int): Flag indicating to have no output-messages from this method
        full (bool): Flag indicating that only fully convertible glycans should be returned, i.e. all modifications
            such as 3-Anhydro-[...] are also present in the SMILES

    Returns:
        Generator generating pairs of type (IUPAC, SMILES) items giving the converted SMILES formulas for the IUPACs.
    """
    logging.basicConfig(level=verbose)

    factory = MonomerFactory()
    glycans = preprocess_glycans(glycan, glycan_list, glycan_file)
    if len(glycans) == 0 and glycan_generator is None:
        logging.info("List of glycans is empty")
        return

    # Convert the glycans ...
    if len(glycans) != 0:
        for glycan in glycans:
            yield generate(glycan, factory, full)

    # Convert the glycans ...
    if glycan_generator is not None:
        for glycan in glycan_generator:
            yield generate(glycan, factory, full)


def generate(glycan, factory, full):
    """
    Actually generate the SMILES string based on the glycan given in IUPAC notation

    Parameters:
        glycan (str): Glycan molecule described by its IUPAC string
        factory (MonomerFactory): Factory to generate the monomers from
        full (bool): flag indicating to only return SMILES string that include all modifications from the IUPAC

    Returns:
        A pair of glycan represented with its IUPAC string and SMILES string
    """
    try:
        # ... by passing them to the glycan class to parse them and return them as intended
        return glycan, Glycan(glycan, factory, full=full).get_smiles()

    # catch any exception at glycan level to not destroy the whole pipeline because of one mis-formed glycan
    except ParseError as e:
        logging.error(f"An exception occurred with {glycan}: {e.__class__}\n"
                      f"Error message: {e.__str__()}")
        return glycan, ""
    except Exception:
        logging.error(f"An unexpected error occurred with with {glycan}. This glycan cannot be parsed.")
        return glycan, ""
