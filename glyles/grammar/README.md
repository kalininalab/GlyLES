# Implemented Grammar in GlyLES

## Implemented Modifications

| Modification     | Example              | parsable | convertable |
|------------------|----------------------|----------|-------------|
| Sulfur groups    | Gal3S                | &check;  | &check;     |
| Phosphate groups | Gal3P                | &check;  | &check;     |
| Fluor groups     | Gal3F                | &check;  | &cross;     |
| Nitrogen gropus  | GalN<br>Gal3N        | &check;  | &check;     |
| Acid groups      | Gal3S<br>2-O-Ac-Gal  | &check;  | &check;     |
| Methyl groups    | Gal2Me<br>2-O-Me-Gal | &check;  | &check;     |
| Desoxygenation   | 3dGal                | &check;  | &cross;     |
| ???              | 3eGal                | &check;  | &cross;     |
| Enantiomers      | L-Gal<br>D-Gal       | &check;  | &cross;     |

###### parsable 
the modification might be in the input glycan and will not raise an exception when parsed. When
converting glycans with parsable modifications, the modification might not necessarily be present in the resulting
SMILES string.

###### convertible
the modification can also be converted into a SMILES string.