grammar Glycan;
start:
    '#' (LBRACK branch RBRACK)* begin '#';
begin:
    branch deriv ' ' TYPE
    | branch deriv
    | deriv ' ' TYPE
    | deriv;
branch:
    deriv con
    | deriv con branch
    | LBRACE branch RBRACE
    | deriv con LBRACE branch RBRACE branch
    | deriv con LBRACE branch RBRACE LBRACE branch RBRACE branch
    | deriv con LBRACE branch RBRACE LBRACE branch RBRACE LBRACE branch RBRACE branch;
deriv:
    modi* saci+ modi* RING? modi* TYPE?;
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
    | DASH END;
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
    '3oxoMyr' | 'A' | 'Ac' | 'Ach' | 'Aep' | 'Al' | 'Ala' | 'Alloc' | 'Allyl' | 'aLnn' | 'Am' | 'Ang' | 'Asp' | 'Beh'
    | 'Bn' | 'Boc' | 'Br' | 'Bu' | 'But' | 'Bz' | 'Cbz' | 'Cct' | 'Cer' | 'Ceroplastic' | 'Cet' | 'Cho' | 'cHx' | 'Cin'
    | 'Cl' | 'ClAc' | 'Cm' | 'Coum' | 'Crt' | 'Cys' | 'DCA' | 'Dce' | 'Dco' | 'DD' | 'Dec' | 'Dhp' | 'DL' | 'DMT'
    | 'Dod' | 'en' | 'eSte' | 'Et' | 'Etg' | 'EtN' | 'Etn' | 'F' | 'Fer' | 'Fmoc' | 'Fo' | 'Gc' | 'Geddic' | 'gLnn'
    | 'Glu' | 'Gly' | 'Gro' | 'Hp' | 'Hpo' | 'Hse' | 'HSer' | 'Hx' | 'Hxo' | 'I' | 'Lac' | 'Lacceroic' | 'Lau' | 'LD'
    | 'Leu' | 'Lev' | 'Lin' | 'LL' | 'Lys' | 'Mal' | 'Mar' | 'Me' | 'Mel' | 'MMT' | 'MOM' | 'Mon' | 'Myr' | 'N3'
    | 'NAP' | 'Ner' | 'Nn' | 'Nno' | 'Non' | 'Ns' | 'Oc' | 'Oco' | 'Ole' | 'oNB' | 'Orn' | 'Pam' | 'Pe' | 'Ph'
    | 'Phthi' | 'Pic' | 'Pico' | 'Piv' | 'PMB' | 'PMP' | 'Poc' | 'Pp' | 'Pr' | 'Pro' | 'Prop' | 'Psyllic' | 'Pyr'
    | 'S' | 'Ser' | 'Sin' | 'Ste' | 'TBDPS' | 'TBS' | 'tBu' | 'TCA' | 'TES' | 'Tf' | 'TFA' | 'THP' | 'Thr' | 'Tig'
    | 'TIPS' | 'TMS' | 'Tr' | 'Troc' | 'Ts' | 'Udo' | 'Ulo' | 'ulo' | 'Und' | 'Vac' | 'Vl';
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
