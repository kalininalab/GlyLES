# GlyLES

Tool to convert IUPAC representation of Glycans into SMILES representation.

## Limitations

This implementation currently only works for glycans that fulfil certain properties:

* The structure of the gylcan is represented as tree of the monomers with maximal branching factor 3.
* Only the 15 implemented monomers may participate in the glycan (see below)

For an overview of the implemented monomers, please look at the README in the grammar folder.

## Workflow

Convert the IUPAC into a SMILES representation using the handy `convert` method

```python
import glyles
glyles.convert(glycan="Man(a1-2)Man", output_file="./blbabl/out.txt")
```

You can also use the `convert` method as a generator
```python
for smiles in glyles.convert(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"], generator=True):\
    print(smiles)
```

In general, the `convert` method supports three different kinds of input.
* single glycan, e.g. `convert(glycan="Man(a1-2)Man)"`,
* a list of glycans, e.g. `convert(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"])`, and
* a file of glycans, e.g. `convert(glycan_file="./glycans.txt")`.<br>Here its important that the file many only 
contain one IUPAC per line.

The output can also be manifold. Supported are 
* `stdout` when specifying no output-related argument, or
* in an `output file`, e.g. `convert(glycan="Man(a1-2)Man", output_file="./out.txt")`.<br>Here each line of the 
output will state the input IUPAC and the output SMILES separated with a comma. Or
* as a `generator`, e.g. `convert(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"], generator=True)`. <br> Then the 
generator will yield the SMILES individually. It is important to note that the `generator` argument overrides the other 
outputs. So, there is always exactly one output method used.</p>
