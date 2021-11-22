# GlyLES

Tool to convert IUPAC representation of Glycans into SMILES representation. This repo is still in the development phase;
so, feel free to report any errors in the issues section.

## Limitations

This implementation currently only works for glycans that fulfil certain properties:

* The structure of the glycan is represented as tree of the monomers with maximal branching factor 3.
* Only the 15 implemented monomers may participate in the glycan (see below)

For an overview of the implemented monomers, please look at the [README](glyles/grammar/README.md) in the grammar folder.

## Workflow

Convert the IUPAC into a SMILES representation using the handy `convert` method

```python
import glyles
glyles.convert(glycan="Man(a1-2)Man", output_file="./blbabl/out.txt")
```

You can also use the `convert_generator` method to get generator over all SMILES:

```python
for smiles in glyles.convert_generator(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"]):
        print(smiles)
```

In general, the `convert` and `convert_generator` methods supports the same for types of input. The samples are shown for `convert` but its the same for `convert_generator`.

* single glycan, e.g. `convert(glycan="Man(a1-2)Man)"`,
* a list of glycans, e.g. `convert(glycan_list=["Man(a1-2)Man a", "Man(a1-2)Man b"])`, and
* a file of glycans, e.g. `convert(glycan_file="./glycans.txt")`.Here its important that the file many only contain one IUPAC per line.
* for better runtime one can also provide a generator as input, e.g. `convert(glycan_generator=some_generator)`

The output also can be manifold for `convert`. For `convert_generator` there is one one output. `convert` supports

* `stdout` when specifying no output-related argument, or
* writing to an `output_file`, e.g. `convert(glycan="Man(a1-2)Man", output_file="./out.txt")`. Here each line of the output will state the input IUPAC and the output SMILES separated with a comma.

In case of `convert_generator` the outputs only contain the SMILES strings in the order of the arguments (first `glycan`, then `glycan_list`, `glycan_file`, and `glycan_generator`).
