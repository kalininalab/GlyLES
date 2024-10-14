grammar Glycan;

start:
    '#' (LBRACK branch RBRACK)* begin '#';
begin:
    branch glycan ' ' TYPE
    | branch glycan
    | glycan ' ' TYPE
    | glycan;
branch:
    glycan con
    | glycan con branch
    | LBRACE branch RBRACE
    | glycan con LBRACE branch RBRACE branch
    | glycan con LBRACE branch RBRACE LBRACE branch RBRACE branch
    | glycan con LBRACE branch RBRACE LBRACE branch RBRACE LBRACE branch RBRACE branch;
glycan:
    deriv RING TYPE | deriv TYPE | deriv RING | deriv;
deriv:
    saci+ RING TYPE modi* | saci+ TYPE modi* | saci+ RING modi* | saci+ modi*
    | modi+ saci+ RING TYPE modi* | modi+ saci+ TYPE modi* | modi+ saci+ RING modi* | modi+ saci+ modi*;
saci: COUNT | SAC;
con:
    LPAR typi NUM DASH qnum RPAR
    | LPAR NUM DASH qnum RPAR
    | typi NUM DASH qnum
    | typi NUM;
add:
    EQ LBRACK ct? NUM (COLON ct? NUM)* RBRACK
    | C LBRACK NUM (COLON NUM)* RBRACK;
fgi:
	COUNT | bridge | FG;
carb:
    (AI | A | I)? CARBON NUM add*;
modi:
    (NUM COLON)? NUM DASH ANHYDRO DASH
    | NUM DASH bridge DASH fgi DASH
    | DASH? NUM? bridge* fgi
    | NUM ((COLON NUM)? D | E | carb)
    | HEAD
    | HEADD DASH
    | DASH END;  // |bridge deleted
qnum:
	QMARK | NUM;
typi:
	TYPE | QMARK;
bridge:
    CARBON | NITROGEN | OXYGEN | PHOSPHOR;
ct:
    C | T;

SAC:
    'Glc' | 'Man' | 'Gal' | 'Gul' | 'Alt' | 'All' | 'Tal' | 'Ido' | 'Qui' | 'Rha' | 'Fuc' | 'Oli' | 'Tyv'
    | 'Abe' | 'Par' | 'Dig' | 'Col' | 'Ara' | 'Lyx' | 'Xyl' | 'Rib' | 'Kdn' | 'Neu' | 'Sia' | 'Pse' | 'Leg' | 'Aci'
    | 'Bac' | 'Kdo' | 'Dha' | 'Mur' | 'Api' | 'Fru' | 'Tag' | 'Sor' | 'Psi' | 'Ery' | 'Thre' | 'Rul' | 'Xul' | 'Unk'
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
    'Ceroplastic' | 'Lacceroic' | '3oxoMyr' | 'Psyllic' | 'Geddic' | 'Allyl' | 'Phthi'
    | 'aLnn' | 'gLnn' | 'eSte' | 'Coum' | 'HSer' | 'Prop'
    | 'Ach' | 'Aep' | 'Ala' | 'Ang' | 'Asp' | 'Beh' | 'But' | 'Cct' | 'Cer' | 'Cet' | 'Cho' | 'Cin' | 'Crt' | 'Cys'
    | 'Dce' | 'Dco' | 'Dec' | 'Dhp' | 'Dod' | 'Etg' | 'EtN' | 'Etn' | 'Fer' | 'Gro' | 'Glu' | 'Gly' | 'Hpo' | 'Hse'
    | 'Hxo' | 'Lac' | 'Lau' | 'Leu' | 'Lin' | 'Lys' | 'Mal' | 'Mar' | 'Mel' | 'Mon' | 'Myr' | 'Ner' | 'Nno' | 'Non'
    | 'Oco' | 'Ole' | 'Orn' | 'Pam' | 'Pro' | 'Pyr' | 'Ser' | 'Sin' | 'Ste' | 'tBu' | 'Thr' | 'Tig' | 'Und' | 'Vac'
    | 'Udo' | 'Ulo' | 'ulo'
    | 'Ac' | 'Am' | 'Bn' | 'Br' | 'Bu' | 'Bz' | 'Cl' | 'Cm' | 'DD' | 'DL' | 'Et' | 'Fo' | 'Gc' | 'Hp' | 'Hx'
    | 'LD' | 'LL' | 'Me' | 'Nn' | 'Oc' | 'Pe' | 'Ph' | 'Pr' | 'Pp' | 'Tf' | 'Tr' | 'Ts' | 'Vl' | 'en'
    | 'A' | 'F' | 'I' | 'S';
ANHYDRO:
    'Anhydro';
HEAD:
    '0d';
HEADD:
    'D' | 'L';
END:
	'ulosaric' | 'ulosonic' | 'uronic' | 'onic' | 'aric' | 'ol';
COUNT:
    'Hep' | 'Hex' | 'Oct' | 'Pen' | 'Suc';
TYPE:
    'a' | 'b';
RING:
    'p' | 'f';
NUM:
    ('1'..'9') ('0'..'9')*;
COLON:
	',';
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
QMARK:
	'?';


// antlr -Dlanguage=Python3 Glycan.g4
