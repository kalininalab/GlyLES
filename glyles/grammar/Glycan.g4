grammar Glycan;

start:
    branch;
branch:
    SAC CON
    | SAC CON branch
    | '[' branch ']'
    | SAC CON '[' branch ']' SAC CON;

SAC:
    'Fuc'
    | 'Gal' | 'Gal3S' | 'Gal6S' | 'GalNAc' | 'GalNAc4S' | 'GalNAc6S'
    | 'Glc' | 'GlcA' | 'GlcNAc' | 'GlcNAc6S'
    | 'Man'
    | 'Neu5Ac' | 'Neu5Gc'
    | 'Tal';
CON:
    '(' TYPE NUM '-' NUM ')';
TYPE:
    'a' | 'b';
NUM: '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9';

// antlr -DLanguage=Python3 Glycan.g4
