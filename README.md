# GlyLES

Tool to convert IUPAC representation of Glycans into SMILES representation

## Limitations

This implementation currently only works for glycans that fulfil certain properties:

* The order of the monomers can be represented as tree with maximal branching factor less than or equal 3.
* Only the 5 implemented monomers may participate in the glycan
* The single monomers may consist of only 1 ring.
* The C1 atom of the root monomer might have the wrong representation of its chirality in SMILES (Suggestion: Leave the
  root monomer unspecified in its chirality)