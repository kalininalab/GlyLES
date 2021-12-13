# GlyLES

A tool to convert IUPAC representation of Glycans into SMILES representation. This repo is still in the development phase;
so, feel free to report any errors in the issues section.

## Installation

So far, this package can only be downloaded from the python package index. So the installation with `pip` is very easy.
Just type

``````shell
pip install glyles
``````

and you're ready to use it as described below.

## Workflow

Convert the IUPAC into a SMILES representation using the handy `convert` method

``````python
from glyles.converter import convert

convert(glycan="Man(a1-2)Man", output_file="./test.txt")
``````

You can also use the `convert_generator` method to get generator over all SMILES:

``````python
from glyles.converter import convert_generator

for smiles in convert_generator(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"]):
    print(smiles)
``````

In general, the `convert` and `convert_generator` methods supports the same for types of input. The samples are shown
for `convert` but it's the same for `convert_generator`.

* single glycan, e.g. `convert(glycan="Man(a1-2)Man)"`,
* a list of glycans, e.g. `convert(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"])`, and
* a file of glycans, e.g. `convert(glycan_file="./glycans.txt")`. Here its important that the file many only contain one
  IUPAC per line.
* for better runtime one can also provide a generator as input, e.g. `convert(glycan_generator=some_generator)`

The output for `convert` can be manifold as well. For `convert_generator` there is one output. `convert` supports

* `stdout` when specifying no output-related argument, or
* writing to an `output_file`, e.g. `convert(glycan="Man(a1-2)Man", output_file="./out.csv")`. Here each line of the
  output will state the input IUPAC and the output SMILES separated with a comma.

In case of `convert_generator` the outputs only contain the SMILES strings in the order of the arguments (first `glycan`
, then `glycan_list`, `glycan_file`, and `glycan_generator`).

## Limitations

This implementation currently only works for glycans that fulfil certain properties:

* The structure of the glycan is represented as tree of the monomers with maximal branching factor 3.
* Only the 23 implemented monomers may participate in the glycan (see below)

For an overview of the implemented monomers, please look at the [README](glyles/grammar/README.md) in the grammar
folder. You can get a python-readable list of the currently included monomers with the following code

`````python
from glyles.glycans.factory import MonomerFactory

MonomerFactory.monomers()
`````

## Poetry

To develop this package I used the poetry package manager (see [here](https://python-poetry.org/) for a detailed
instruction). It has basically the same functionality as conda but supports the package management i/m.m.o. better and
also gives to opportunity to distinguish packages into those that are needed to use the package and those that are
needed in the development of the package. In order to enable others to work on this repository, I also published the
exact specifications of my poetry environment.
