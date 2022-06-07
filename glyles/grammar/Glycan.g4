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
    MOD* SAC+ RING TYPE MOD* | MOD* SAC+ TYPE MOD* | MOD* SAC+ RING MOD* | MOD* SAC+ MOD*;

SAC:
    COUNT | 'Glc' | 'Man' | 'Gal' | 'Gul' | 'Alt' | 'All' | 'Tal' | 'Ido' | 'Qui' | 'Rha' | 'Fuc' | 'Oli'
    | 'Tyv' | 'Abe' | 'Par' | 'Dig' | 'Col' | 'Ara' | 'Lyx' | 'Xyl' | 'Rib' | 'Kdn' | 'Neu' | 'Pse' | 'Leg' | 'Aci'
    | 'Bac' | 'LDManHep' | 'Kdo' | 'Dha' | 'DDManHep' | 'Mur' | 'Api' | 'Fru' | 'Tag' | 'Sor' | 'Psi' | 'Ery';
MOD:
    NUM ',' NUM '-Anhydro-'
    | NUM '-O-' FG '-'
    | NUM ('N' | 'O' | 'C') FG
    | (NUM | 'O' | 'N' | 'C') FG
    | NUM ('d' | 'e')
    | '-ulosaric' | '-ulosonic' | '-uronic' | '-onic' | '-aric' | '-ol'
    | '0d' | 'D-' | 'L-' | FG;
FG:
    COUNT | 'SerAc' | 'ThrAc'
    | 'Coum' | 'Prop'
    | 'Ach' | 'Aep' | 'Ala' | 'Asp' | 'Beh' | 'But' | 'Cho' | 'Cin' | 'Cys' | 'Dco' | 'Dhp' | 'Etg' | 'Etn' | 'Fer'
    | 'Gro' | 'Glu' | 'Gly' | 'Hxo' | 'Lac' | 'Lau' | 'Leu' | 'Lin' | 'Lys' | 'Mal' | 'Mar' | 'Myr' | 'Oco' | 'Ole'
    | 'Orn' | 'Pam' | 'Pro' | 'Pyr' | 'Ser' | 'Sin' | 'Ste' | 'Suc' | 'Thr' | 'Vac' | 'Ulo' | 'ulo'
    | 'Ac' | 'Am' | 'Bn' | 'Br' | 'Bz' | 'Cl' | 'Cm' | 'Et' | 'Fo' | 'Gc' | 'Gr' | 'Me' | 'Oc' | 'Ph' | 'Pp' | 'Tf'
    | 'Tr' | 'Ts' | 'Vl' | 'en'
    | 'A' | 'N' | 'F' | 'I' | 'S' | 'P';
COUNT:
    'Hep' | 'Hex' | 'Oct' | 'Pen';
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
