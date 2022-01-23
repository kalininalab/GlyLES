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
    deriv RING TYPE | deriv TYPE | deriv RING | deriv;
deriv:
    MOD* SAC MOD*;

SAC:
    'Glc' | 'Man' | 'Gal' | 'Gul' | 'Alt' | 'All' | 'Tal' | 'Ido' | 'Qui' | 'Rha' | 'Fuc' | 'Oli' | 'Tyv' | 'Abe'
    | 'Par' | 'Dig' | 'Col' | 'Ara' | 'Lyx' | 'Xyl' | 'Rib' | 'Kdn' | 'Neu' | 'Pse' | 'Leg' | 'Aci' | 'Bac'
    | 'LDmanHep' | 'Kdo' | 'Dha' | 'DDmanHep' | 'Mur' | 'Api' | 'Fru' | 'Tag' | 'Sor' | 'Psi';
MOD:
    NUM '-O-Me-' | NUM '-O-Ac-' | 'NAc' | NUM 'Me' | NUM 'Ac' | NUM 'd' | NUM 'S' | NUM 'P' | NUM 'N' | 'Ac' | 'A' | 'N';
CON:
    '(' TYPE NUM '-' NUM ')';
TYPE:
    'a' | 'b';
RING:
    'p' | 'f';
NUM:
    '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9';  // [1-9] ??

// antlr -DLanguage=Python3 Glycan.g4
