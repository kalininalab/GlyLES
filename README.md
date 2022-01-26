# GlyLES

A tool to convert IUPAC representation of Glycans into SMILES representation. This repo is still in the development 
phase; so, feel free to report any errors or issues.

## Specification and (current) Limitations

The exact specification, we're referring to when talking about "IUPAC representations of glycan", is given in the 
"Notes" section of this [website](https://www.ncbi.nlm.nih.gov/glycans/snfg.html). But as this package is still in the 
development phase, not everything of the specification is implemented yet (especially not all monomers and side chains 
you can attach to monomers).

This implementation currently only works for glycans that fulfil certain properties:

* Linkages have to be explicit, i.e. `(a1-4)`
* The structure of the glycan can be represented as tree of the monomers with maximal branching factor 2.
* All root monomers (e.g. Glc, but not GlcNAc) from this [website](https://www.ncbi.nlm.nih.gov/glycans/snfg.html)
* Some modifications can be added to the monomers, please see the [README](glyles/grammar/README.md) in the grammar
folder for more information on this. 

## Installation

So far, this package can only be downloaded from the python package index. So the installation with `pip` is very easy.
Just type

``````shell
pip install glyles
``````

and you're ready to use it as described below. Use 

````shell
pip install --upgrade glyles
````

to upgrade the glyles package to the most recent version.

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
* a list of glycans, e.g. `convert(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man"])`, and
* a file of glycans, e.g. `convert(glycan_file="./glycans.txt")`. Here its important that the file many only contain one
  IUPAC per line.
* for better runtime one can also provide a generator as input, e.g. `convert(glycan_generator=some_generator)`

The output for `convert` can be manifold as well. For `convert_generator` there is one output. `convert` supports

* `stdout` when specifying no output-related argument, or
* returned as list of tuples if `returning=true` is set, or
* writing to an `output_file`, e.g. `convert(glycan="Man(a1-2)Man", output_file="./out.csv")`.

Here each line of the output will state the input IUPAC and the output SMILES. In case of written output separated with 
a comma.


## Poetry

To develop this package I used the poetry package manager (see [here](https://python-poetry.org/) for a detailed
instruction). It has basically the same functionality as conda but supports the package management better and also 
supports distinguishing packages into those that are needed to use the package and those that are needed in the 
development of the package. In order to enable others to work on this repository, we also publish the exact 
specifications of out poetry environment.
