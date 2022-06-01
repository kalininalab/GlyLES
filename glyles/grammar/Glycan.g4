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
    | glycan CON '[' branch ']' branch
    | glycan CON '[' branch ']' '[' branch ']' branch;
glycan:
    deriv RING TYPE | deriv TYPE | deriv RING | deriv;
deriv:
    MOD* SAC RING TYPE MOD* | MOD* SAC TYPE MOD* | MOD* SAC RING MOD* | MOD* SAC MOD*;

SAC:
    'Glc' | 'Man' | 'Gal' | 'Gul' | 'Alt' | 'All' | 'Tal' | 'Ido' | 'Qui' | 'Rha' | 'Fuc' | 'Oli'
    | 'Tyv' | 'Abe' | 'Par' | 'Dig' | 'Col' | 'Ara' | 'Lyx' | 'Xyl' | 'Rib' | 'Kdn' | 'Neu' | 'Pse' | 'Leg' | 'Aci'
    | 'Bac' | 'LDManHep' | 'Kdo' | 'Dha' | 'DDManHep' | 'Mur' | 'Api' | 'Fru' | 'Tag' | 'Sor' | 'Psi' | 'Ery';
MOD:
    NUM ',' NUM '-Anhydro-' | NUM '-O-' ('Ac' | 'Bn' | 'Bz' | 'Et' | 'Me') '-'
    | NUM 'N' ('SerAc' | 'ThrAc' | 'But' | 'Me')
    | (NUM | 'O') ('Prop'
        | 'Ala' | 'Beh' | 'But' | 'Cho' | 'Etn' | 'Fer' | 'Lac' | 'Gro' | 'Hxo' | 'Lau' | 'Lin' | 'Oco' | 'Ole' | 'Pyr'
        | 'Ac' | 'Am' | 'Bn' | 'Bz' | 'Et' | 'Fo' | 'Me' | 'Gc' | 'Gr' | 'Ph' | 'Tf' | 'Tr'
        | 'Ts' | 'S' | 'P' | 'N' | 'F' | 'd' | 'e')
    | 'N' ('But' | 'Lac' | 'Ac' | 'Bz' | 'Fo' | 'Gc' | 'Me')
    | '-ulosaric' | '-uronic' | '-onic' | '-aric' | '-ol'
    | 'Ala' | 'Hep' | 'Hex' | 'Ser' | 'Pen'
    | '0d' | 'D-' | 'L-' | 'Ac' | 'CN'
    | 'A' | 'N' | 'P';
CON:
    '(' TYPE NUM '-' NUM ')'
    | TYPE NUM '-' NUM
    | TYPE NUM;
TYPE:
    'a' | 'b';
RING:
    'p' | 'f';
NUM:
    '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9';  // [1-9] ??

// antlr -DLanguage=Python3 Glycan.g4
