# GlyLES
Tool to convert IUPAC representation of Glycans into SMILES representation

## Important
This implementation currently only works for glycans that fulfil certain properties:
 * The order of the monomers can be represented as tree with maximal branching factor less than or equal 2.
 * Only the 5 implemented monomers may participate in the glycan
 * The single monomers may consist of only 1 ring.