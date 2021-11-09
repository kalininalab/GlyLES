grammar Glycan;

start:
    branch;
branch:
    SAC CON
    | SAC CON branch
    | '[' branch ']'
    | SAC CON '[' branch ']' SAC CON;

CON:
    '(' TYPE NUM '-' NUM ')';
SAC:
    'Glc' | 'Glu' | 'Tal' | 'Man' | 'Gal';
TYPE:
    'a' | 'b';
NUM: '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9';

// antlr -DLanguage=Python3 Glycan.g4
