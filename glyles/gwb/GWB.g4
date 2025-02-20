grammar GWB;

start: 
    bos branch (RBRACK branch)? DOLLAR;
branch:
    con deriv
    | LPAR branch RPAR branch
    | con deriv (LPAR branch RPAR)? branch;
deriv:
    (HEADD DASH)? modi* saci+ modi* (COMMA RING)? (COMMA? AT NUM)? 
    | modi* (COMMA? AT NUM)?;
saci: 
    COUNT | SAC;
con: 
    DOUBLEDASH (qnum typi)? qnum?;
add:
    EQ LBRACK ct? NUM (COMMA ct? NUM)* RBRACK
    | C LBRACK NUM (COMMA NUM)* RBRACK;
fgi:
	COUNT | bridge | FG;
carb:
    (AI | A | I)? CARBON NUM add*;
modi:
    (NUM COMMA)? NUM DASH ANHYDRO DASH
    | NUM DASH bridge DASH fgi DASH
    | DASH? NUM? bridge* fgi
    | QMARK NUM FG
    | NUM ((COMMA NUM)? D | E | U | carb)
    | HEAD
    | HEADD DASH
    | DASH END;
qnum:
	QMARK | NUM (SLASH NUM)*;
typi: 
    TYPE | QMARK;
bridge:
    CARBON | NITROGEN | OXYGEN | PHOSPHOR;
ct:
    C | T;
bos: 
    FREEEND | REDEND;







// Tokens
FREEEND:
    'freeEnd';
REDEND:
    'redEnd';
SAC:
    'Glc' | 'Man' | 'Gal' | 'Gul' | 'Alt' | 'All' | 'Tal' | 'Ido' | 'Qui' | 'Rha' | 'Fuc' | 'Oli' | 'Tyv'
    | 'Abe' | 'Par' | 'Dig' | 'Col' | 'Ara' | 'Lyx' | 'Xyl' | 'Rib' | 'Kdn' | 'KDN' | 'Neu' | 'Sia' | 'Pse' | 'Leg' | 'Aci'
    | 'Bac' | 'Kdo' | 'KDO' | 'Dha' | 'Mur' | 'Api' | 'Fru' | 'Tag' | 'Sor' | 'Psi' | 'Ery' | 'Thre' | 'Rul' | 'Xul' | 'Unk'
    | 'Ace' | 'Aco' | 'Asc' | 'Fus' | 'Ins' | 'Ko' | 'Pau' | 'Per'| 'Sed' | 'Sug' | 'Vio' | 'Xlu' | 'Yer' | 'Erwiniose';
CARBON:
	'C';
NITROGEN:
    'N';
OXYGEN:
    'O';
PHOSPHOR:
    'P';
FG:
    'Ceroplastic' | 'Lacceroic' | '3oxoMyr' | 'Psyllic' | 'Geddic' | 'Alloc' | 'Allyl' | 'Phthi' | 'TBDPS' | 'aLnn'
    | 'ClAc' | 'Coum' | 'eSte' | 'Fmoc' | 'gLnn' | 'HSer' | 'Pico' | 'Prop' | 'TIPS' | 'triN' | 'Troc' | 'Ach' | 'Aep'
    | 'Ala' | 'Ang' | 'Asp' | 'Beh' | 'Boc' | 'But' | 'Cbz' | 'Cct' | 'Cer' | 'Cet' | 'Cho' | 'cHx' | 'Cin' | 'Crt'
    | 'Cys' | 'DCA' | 'Dce' | 'Dco' | 'Dec' | 'Dhp' | 'DMT' | 'Dod' | 'Etg' | 'Etn' | 'EtN' | 'Fer' | 'Glu' | 'Gly'
    | 'Gro' | 'Hpo' | 'Hse' | 'Hxo' | 'Lac' | 'Lau' | 'Leu' | 'Lev' | 'Lin' | 'Lys' | 'Mal' | 'Mar' | 'Mel' | 'MMT'
    | 'MOM' | 'Mon' | 'Myr' | 'NAP' | 'Ner' | 'Nno' | 'Non' | 'Oco' | 'Ole' | 'oNB' | 'Orn' | 'Pam' | 'Pic' | 'Piv'
    | 'PMB' | 'PMP' | 'Poc' | 'Pro' | 'Pyr' | 'Ser' | 'Sin' | 'Ste' | 'TBS' | 'tBu' | 'TCA' | 'TES' | 'TFA' | 'THP'
    | 'Thr' | 'Tig' | 'TMS' | 'Udo' | 'Ulo' | 'ulo' | 'Und' | 'Vac' | 'Ac' | 'Al' | 'Am' | 'Bn' | 'Br' | 'Bu' | 'Bz'
    | 'Cl' | 'Cm' | 'DD' | 'DL' | 'en' | 'Et' | 'Fo' | 'Gc' | 'Hp' | 'Hx' | 'LD' | 'LL' | 'Me' | 'Nn' | 'Ns'
    | 'Oc' | 'Pe' | 'Ph' | 'Pp' | 'Pr' | 'Tf' | 'Tr' | 'Ts' | 'Vl' | 'A' | 'F' | 'I' | 'S';
ANHYDRO:
    'Anhydro';
HEAD:
    '0d';
HEADD:
    'D' | 'L';
END:
	'ulosaric' | 'ulosonic' | 'uronic' | 'onic' | 'aric' | 'ol';
COUNT:
    'Hept' | 'Hep' | 'Hex' | 'Oct' | 'Pen' | 'Suc';
TYPE: 
    'a' | 'b';
RING:
    'p' | 'f' | 'o';
NUM:
    ('1'..'9') ('0'..'9')*;
AT:
    '@';
COMMA:
	',';
SEMICOLON: 
    ';';
DOUBLEDASH: 
    '--';
DASH: 
    '-';
LPAR: 
    '(';
RPAR: 
    ')';
LBRACE:
	'[';
RBRACE:
	']';
LBRACK:
    '{';
RBRACK:
    '}';
A:
	'a';
AI:
	'ai';
C:
	'c';
D:
	'd';
E:
	'e';
T:
    't';
EQ:
	'=';
I:
	'i';
U:
    'u';
SLASH:
    '/';
DOLLAR:
    '$';
HASH:
    '#';
SPACE:
    ' ';
QMARK: 
    '?';

// antlr -Dlanguage=Python3 GWB.g4
