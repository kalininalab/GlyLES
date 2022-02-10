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
    MOD* SAC RING TYPE MOD* | MOD* SAC TYPE MOD* | MOD* SAC RING MOD* | MOD* SAC MOD*;

SAC:
    'Hex' | 'Pen' | 'Glc' | 'Man' | 'Gal' | 'Gul' | 'Alt' | 'All' | 'Tal' | 'Ido' | 'Qui' | 'Rha' | 'Fuc' | 'Oli'
    | 'Tyv' | 'Abe' | 'Par' | 'Dig' | 'Col' | 'Ara' | 'Lyx' | 'Xyl' | 'Rib' | 'Kdn' | 'Neu' | 'Pse' | 'Leg' | 'Aci'
    | 'Bac' | 'LDmanHep' | 'Kdo' | 'Dha' | 'DDmanHep' | 'Mur' | 'Api' | 'Fru' | 'Tag' | 'Sor' | 'Psi';
MOD:
    NUM ',' NUM '-Anhydro-' | NUM ('-O-Ac-' | '-O-Me-')
    | NUM 'N' ('SerAc' | 'ThrAc' | 'But' | 'Me')
    | (NUM | 'O') ('Prop' | 'Ala' | 'But'| 'Etn' | 'Fer' | 'Lac' | 'Gro' | 'Hxo' | 'Lin' | 'Oco' | 'Ole' | 'Pyr'
        | 'Ac' | 'Am' | 'Fo' | 'Me' | 'Gc' | 'Gr' | 'PP' | 'S' | 'P' | 'N' | 'F' | 'd' | 'e')
    | 'N' ('But' | 'Lac' | 'Ac' | 'Fo' | 'Me')
    | '-ulosaric' | '-uronic' | '-onic' | '-aric' | 'Ala' | 'Ser' | '-ol' | '0d' | 'D-' | 'L-' | 'Ac' | 'CN' | 'A' | 'N';
CON:
    '(' TYPE NUM '-' NUM ')';
TYPE:
    'a' | 'b';
RING:
    'p' | 'f';
NUM:
    '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9';  // [1-9] ??

// antlr -DLanguage=Python3 Glycan.g4
