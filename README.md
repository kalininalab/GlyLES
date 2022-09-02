# GlyLES

![testing](https://github.com/kalininalab/glyles/actions/workflows/test.yaml/badge.svg) ![piwheels](https://img.shields.io/piwheels/v/glyles) ![PyPI - Downloads](https://img.shields.io/pypi/dm/glyles) [![codecov](https://codecov.io/gh/kalininalab/GlyLES/branch/main/graph/badge.svg)](https://codecov.io/gh/kalininalab/glyles)

A tool to convert IUPAC representation of Glycans into SMILES representation. This repo is still in the development 
phase; so, feel free to report any errors or issues.

## Specification and (current) Limitations

The exact specification we're referring to when talking about "IUPAC representations of glycan", is given in the 
"Notes" section of this [website](https://www.ncbi.nlm.nih.gov/glycans/snfg.html). But as this package is still in the 
development phase, not everything of the specification is implemented yet (especially not all side chains you can 
attach to monomers).

This implementation currently only works for glycans that fulfill certain properties:

* The structure of the glycan can be represented as a tree of the monosaccharides with maximal branching factor 3.
* Only basic monomers (e.g. ``Glc``, but not ``GlcNAc``) from this [website](https://www.ncbi.nlm.nih.gov/glycans/snfg.html) 
  (GalNAc is seen as modification of galactose). In addition to this list, Erythrose is parsable.
* Some modifications can be added to the monomers, please see the [README](glyles/grammar/README.md) in the grammar
folder for more information on this. <br>E.g. ``GlcNAc`` is seen to be a modified version of ``Glc``. 

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

## Usage

Convert the IUPAC into a SMILES representation using the handy `convert` method

``````python
from glyles.converter import convert

convert(glycan="Man(a1-2)Man", output_file="./test.txt")
``````

You can also use the `convert_generator` method to get a generator for all SMILES:

``````python
from glyles.converter import convert_generator

for smiles in convert_generator(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"]):
    print(smiles)
``````

In general, the `convert` and `convert_generator` methods support the same types of input. The samples are shown
for `convert` but it's the same for `convert_generator`.

* single glycan, e.g. `convert(glycan="Man(a1-2)Man)"`,
* a list of glycans, e.g. `convert(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man"])`, and
* a file of glycans, e.g. `convert(glycan_file="./glycans.txt")`. Here its important that the file many only contain one
  IUPAC per line.
* for better runtime one can also provide a generator as input, e.g. `convert(glycan_generator=some_generator)`

Any output consists of tuples of the form `(input_iupac, smiles)`. `convert_generator` returns those tuples in a 
generator. `convert`-output is a bit different:

* By default, the method returns the list of tuples.
* You can print them to `stdout` by specifying `returning=False`.
* The tuples can also be written `output_file`, e.g. `convert(glycan="Man(a1-2)Man", output_file="./out.csv")`.

## Notation of glycans

There are multiple different notations for glycans in IUPAC. So, according to the 
[SNGF specification](https://www.ncbi.nlm.nih.gov/glycans/snfg.html), `Man(a1-4)Gal`, `Mana1-4Gal`, and `Mana4Gal` 
all describe the same disaccharide. This is also covered in this package as all three notations will be parsed into the 
same tree of monosaccharides and result in the same SMILES string.


## Poetry

To develop this package, we use the poetry package manager (see [here](https://python-poetry.org/) for detailed
instruction). It has basically the same functionality as conda but supports the package management better and also 
supports distinguishing packages into those that are needed to use the package and those that are needed in the 
development of the package. To enable others to work on this repository, we also publish the exact 
specifications of our poetry environment.
