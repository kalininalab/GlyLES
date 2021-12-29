grammar Glycan;

start:
    '{' branch glycan ' ' TYPE '}'
    | '{' branch glycan '}'
    | '{' glycan ' ' TYPE '}'
    | '{' glycan '}';
branch:
    glycan CON
    | glycan CON branch
    | '[' branch ']'
    | glycan CON '[' branch ']' branch;

glycan:
    SAC RING TYPE | SAC TYPE | SAC RING | SAC;

SAC:
    'Glc' | 'Glc4S' | 'Glc6S' | 'GlcA' | 'GlcN' | 'GlcNAc' | 'GlcNAc3S' | 'GlcNAc6S'
    | 'Man'
    | 'Gal' | 'Gal3S' | 'Gal3S4S' | 'Gal3S6S' | 'Gal4S' | 'Gal4S6S' | 'Gal6S' | 'GalNAc' | 'GalNAc3S' | 'GalNAc4S'
    | 'GalNAc4S6S' | 'GalNAc6S'
    | 'Gul' | 'Alt' | 'All' | 'Tal' | 'Ido' | 'Qui' | 'Rha' | 'Fuc' | 'Oli' | 'Tyv' | 'Abe' | 'Par' | 'Dig' | 'Col'
    | 'Ara' | 'Lyx' | 'Xyl' | 'Rib' | 'Kdn' | 'Neu' | 'Pse' | 'Leg' | 'Aci' | 'Bac' | 'LDmanHep' | 'Kdo' | 'Dha'
    | 'DDmanHep' | 'Mur' | 'Api' | 'Fru' | 'Tag' | 'Sor' | 'Psi';
CON:
    '(' TYPE NUM '-' NUM ')';
TYPE:
    'a' | 'b';
RING:
    'p' | 'f';
NUM:
    '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9';  // [1-9] ??

// antlr -DLanguage=Python3 Glycan.g4
