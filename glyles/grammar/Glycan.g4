grammar Glycan;

start:
    branch SAC ' ' TYPE
    | branch SAC
    | SAC ' ' TYPE
    | SAC;
branch:
    SAC CON
    | SAC CON branch
    | '[' branch ']'
    | SAC CON '[' branch ']' branch;

SAC:
    'Fuc'
    | 'Gal' | 'Gal3S' | 'Gal3S4S' | 'Gal3S6S' | 'Gal4S' | 'Gal4S6S' | 'Gal6S' | 'GalNAc' | 'GalNAc3S' | 'GalNAc4S'
    | 'GalNAc4S6S' | 'GalNAc6S'
    | 'Glc' | 'Glc4S' | 'Glc6S' | 'GlcA' | 'GlcN' | 'GlcNAc' | 'GlcNAc3S' | 'GlcNAc6S'
//    | 'Kdn'
    | 'Man'
//    | 'Neu5Ac' | 'Neu5Gc'
    | 'Tal';
CON:
    '(' TYPE NUM '-' NUM ')';
TYPE:
    'a' | 'b';
NUM:
    '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9';  // [1-9] ??

// antlr -DLanguage=Python3 Glycan.g4
