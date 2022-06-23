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
    | glycan CON '[' branch ']' '[' branch ']' branch
    | glycan CON '[' branch ']' '[' branch ']' '[' branch ']' branch;
glycan:
    deriv RING TYPE | deriv TYPE | deriv RING | deriv;
deriv:
    SAC+ RING TYPE MOD* | SAC+ TYPE MOD* | SAC+ RING MOD* | SAC+ MOD*
    | MOD+ SAC+ RING TYPE MOD* | MOD+ SAC+ TYPE MOD* | MOD+ SAC+ RING MOD* | MOD+ SAC+ MOD*;

SAC:
    COUNT | 'Glc' | 'Man' | 'Gal' | 'Gul' | 'Alt' | 'All' | 'Tal' | 'Ido' | 'Qui' | 'Rha' | 'Fuc' | 'Oli' | 'Tyv'
    | 'Abe' | 'Par' | 'Dig' | 'Col' | 'Ara' | 'Lyx' | 'Xyl' | 'Rib' | 'Kdn' | 'Neu' | 'Sia' | 'Pse' | 'Leg' | 'Aci'
    | 'Bac' | 'Kdo' | 'Dha' | 'Mur' | 'Api' | 'Fru' | 'Tag' | 'Sor' | 'Psi' | 'Ery' | 'Thre' | 'Rul' | 'Xul'
    | 'Ace' | 'Aco' | 'Asc' | 'Fus' | 'Ins' | 'Ko' | 'Pau' | 'Per'| 'Sed' | 'Sug' | 'Vio' | 'Xlu' | 'Yer' | 'Erwiniose';
MOD:
    NUM ',' NUM '-Anhydro-' | NUM '-Anhydro-'
    | NUM '-O-' FG '-'
    | NUM ('NH' | 'N' | 'O' | 'C') FG
    | (NUM | 'O' | 'NH' | 'N' | 'C') FG
    | NUM ('d' | 'e') | NUM ',' NUM 'd'
    | '-ulosaric' | '-ulosonic' | '-uronic' | '-onic' | '-aric' | '-ol'
    | '0d' | 'D-' | 'L-' | FG;
FG:
    COUNT | 'Coum' | 'Prop'
    | 'Ach' | 'Aep' | 'Ala' | 'Asp' | 'Beh' | 'But' | 'Cho' | 'Cer' | 'Cet' | 'Cin' | 'Cys' | 'Dco' | 'Dhp' | 'Etg'
    | 'Etn' | 'Fer' | 'Gro' | 'Glu' | 'Gly' | 'Hse' | 'Hxo' | 'Lac' | 'Lau' | 'Leu' | 'Lin' | 'Lys' | 'Mal' | 'Mar'
    | 'Myr' | 'Non' | 'Oco' | 'Ole' | 'Orn' | 'Pam' | 'Pro' | 'Pyr' | 'Ser' | 'Sin' | 'Ste' | 'Thr' | 'Vac' | 'Ulo'
    | 'ulo'
    | 'Ac' | 'Am' | 'Bn' | 'Br' | 'Bz' | 'Cl' | 'Cm' | 'DD' | 'DL' | 'DL-' | 'Et' | 'Fo' | 'Gc' | 'LD' | 'LL' | 'Me'
    | 'Oc' | 'Ph' | 'Pp' | 'Tf' | 'Tr' | 'Ts' | 'Vl' | 'en'
    | 'A' | 'N' | 'F' | 'I' | 'S' | 'P';
COUNT:
    'Hep' | 'Hex' | 'Oct' | 'Pen' | 'Suc';
CON:
    '(' TYPE NUM '-' NUM ')'
    | '(' NUM '-' NUM ')'
    | TYPE NUM '-' NUM
    | TYPE NUM;
    // | '(' 'z' NUM '-' 'z' ')' | 'z' NUM '-' 'z' | 'z' NUM;
TYPE:
    'a' | 'b' | 'z';
RING:
    'p' | 'f';
NUM:
    '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | '10' | '11' | '12' | '13' | 'z';  // [1-9] ??

// antlr -DLanguage=Python3 Glycan.g4
