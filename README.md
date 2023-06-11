# GlyLES

![testing](https://github.com/kalininalab/glyles/actions/workflows/test.yaml/badge.svg)
[![docs-image](https://readthedocs.org/projects/glyles/badge/?version=latest)](https://glyles.readthedocs.io/en/latest/)
[![piwheels](https://img.shields.io/piwheels/v/glyles)](https://pypi.org/project/glyles/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/glyles)](https://pypi.org/project/glyles/)
[![codecov](https://codecov.io/gh/kalininalab/GlyLES/branch/main/graph/badge.svg)](https://codecov.io/gh/kalininalab/glyles)
[![DOI](https://zenodo.org/badge/431874597.svg)](https://zenodo.org/badge/latestdoi/431874597)

A tool to convert IUPAC representation of Glycans into SMILES representation. This repo is still in the development 
phase; so, feel free to report any errors or issues. The code is available on 
[github](https://github.com/kalininalab/GlyLES/) and the documentation can be found on 
[ReadTheDocs](https://glyles.readthedocs.io/en/latest/index.html).

## Specification and (current) Limitations

The exact specification we're referring to when talking about "IUPAC representations of glycan" or "IUPAC-condensed", 
is given in the "Notes" section of this [website](https://www.ncbi.nlm.nih.gov/glycans/snfg.html). But as this package 
is still in the development phase, not everything of the specification is implemented yet (especially not all side 
chains you can attach to monomers). The structure of the glycan can be represented as a tree of the monosaccharides 
with maximal branching factor 4, i.e., each monomer in the glycan has at most 4 children.

## Installation

So far, this package can only be downloaded from the python package index. So the installation with `pip` is very easy.
Just type

``````shell
pip install glyles
``````

and you're ready to use it as described below. Use 

``````shell
pip install --upgrade glyles
``````

to upgrade the glyles package to the most recent version.

## Basic Usage

### As a Python Package

Convert the IUPAC into a SMILES representation using the handy `convert` method

``````python
from glyles import convert

convert(glycan="Man(a1-2)Man", output_file="./test.txt")
``````

You can also use the `convert_generator` method to get a generator for all SMILES:

``````python
from glyles import convert_generator

for smiles in convert_generator(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"]):
    print(smiles)
``````

For more examples of how to use this package, please see the notebooks in the 
[examples](https://github.com/kalininalab/GlyLES/tree/dev/examples) folder and checkout the documentation on 
[ReadTheDocs](https://glyles.readthedocs.io/en/latest/index.html).

### In the Commandline

As of version 0.5.9, there is a commandline interface to GlyLES which is automatically installed when installing GlyLES 
through pip. The CLI is open for one or multiple IUPAC inputs as individual arguments. Due to the syntax of the 
IUPAC-condensed notation and the argument parsing in commandlines, the IUPAC strings must be given in quotes.

``````shell
glyles -i "Man(a1-2)Man" -o test_output.txt
glyles -i "Man(a1-2)Man" "Fuc(a1-6)Glc" -o test_output.txt
``````

File-input is also possible.
``````shell
glyles -i input_file.txt -o test_output.txt
``````

Providing multiple files and IUPAC-condensed names is als supported.
``````shell
glyles -i input_file1.txt "Man(a1-2)Man" input_file2.txt input_file13.txt "Fuc(a1-6)Glc" -o test_output.txt
``````

## Notation of glycans

There are multiple different notations for glycans in IUPAC. So, according to the 
[SNGF specification](https://www.ncbi.nlm.nih.gov/glycans/snfg.html), `Man(a1-4)Gal`, `Mana1-4Gal`, and `Mana4Gal` 
all describe the same disaccharide. This is also covered in this package as all three notations will be parsed into the 
same tree of monosaccharides and result in the same SMILES string.

This is also described more detailed in a section on [ReadTheDocs](https://glyles.readthedocs.io/en/latest/notes/notation.html).

## Poetry

To develop this package, we use the poetry package manager (see [here](https://python-poetry.org/) for detailed
instruction). It has basically the same functionality as conda but supports the package management better and also 
supports distinguishing packages into those that are needed to use the package and those that are needed in the 
development of the package. To enable others to work on this repository, we also publish the exact 
specifications of our poetry environment.

## Citation

If you use GlyLES in your work, please cite
```
@article{joeres2023glyles,
  title={GlyLES: Grammar-based Parsing of Glycans from IUPAC-condensed to SMILES},
  author={Joeres, Roman and Bojar, Daniel and Kalinina, Olga V},
  journal={Journal of Cheminformatics},
  volume={15},
  number={1},
  pages={1--11},
  year={2023},
  publisher={BioMed Central}
}
```
