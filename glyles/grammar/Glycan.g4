grammar Glycan;

start:
    branch;
branch:
    SAC CON
    | SAC CON branch
    | '[' branch ']'
    | SAC CON '[' branch ']' SAC CON;

SAC:
    'Glc' | 'Glu' | 'Fru' | 'Man' | 'Gal';
CON:
    '(a1-1)' | '(a1-2)' | '(a1-3)' | '(a1-4)' | '(a1-5)' | '(a1-6)' | '(a1-7)' | '(a1-8)' | '(a1-9)'
	| '(a2-1)' | '(a2-2)' | '(a2-3)' | '(a2-4)' | '(a2-5)' | '(a2-6)' | '(a2-7)' | '(a2-8)' | '(a2-9)'
	| '(a3-1)' | '(a3-2)' | '(a3-3)' | '(a3-4)' | '(a3-5)' | '(a3-6)' | '(a3-7)' | '(a3-8)' | '(a3-9)'
	| '(a4-1)' | '(a4-2)' | '(a4-3)' | '(a4-4)' | '(a4-5)' | '(a4-6)' | '(a4-7)' | '(a4-8)' | '(a4-9)'
	| '(a5-1)' | '(a5-2)' | '(a5-3)' | '(a5-4)' | '(a5-5)' | '(a5-6)' | '(a5-7)' | '(a5-8)' | '(a5-9)'
	| '(a6-1)' | '(a6-2)' | '(a6-3)' | '(a6-4)' | '(a6-5)' | '(a6-6)' | '(a6-7)' | '(a6-8)' | '(a6-9)'
	| '(a7-1)' | '(a7-2)' | '(a7-3)' | '(a7-4)' | '(a7-5)' | '(a7-6)' | '(a7-7)' | '(a7-8)' | '(a7-9)'
	| '(a8-1)' | '(a8-2)' | '(a8-3)' | '(a8-4)' | '(a8-5)' | '(a8-6)' | '(a8-7)' | '(a8-8)' | '(a8-9)'
	| '(a9-1)' | '(a9-2)' | '(a9-3)' | '(a9-4)' | '(a9-5)' | '(a9-6)' | '(a9-7)' | '(a9-8)' | '(a9-9)'
	| '(b1-1)' | '(b1-2)' | '(b1-3)' | '(b1-4)' | '(b1-5)' | '(b1-6)' | '(b1-7)' | '(b1-8)' | '(b1-9)'
	| '(b2-1)' | '(b2-2)' | '(b2-3)' | '(b2-4)' | '(b2-5)' | '(b2-6)' | '(b2-7)' | '(b2-8)' | '(b2-9)'
	| '(b3-1)' | '(b3-2)' | '(b3-3)' | '(b3-4)' | '(b3-5)' | '(b3-6)' | '(b3-7)' | '(b3-8)' | '(b3-9)'
	| '(b4-1)' | '(b4-2)' | '(b4-3)' | '(b4-4)' | '(b4-5)' | '(b4-6)' | '(b4-7)' | '(b4-8)' | '(b4-9)'
	| '(b5-1)' | '(b5-2)' | '(b5-3)' | '(b5-4)' | '(b5-5)' | '(b5-6)' | '(b5-7)' | '(b5-8)' | '(b5-9)'
	| '(b6-1)' | '(b6-2)' | '(b6-3)' | '(b6-4)' | '(b6-5)' | '(b6-6)' | '(b6-7)' | '(b6-8)' | '(b6-9)'
	| '(b7-1)' | '(b7-2)' | '(b7-3)' | '(b7-4)' | '(b7-5)' | '(b7-6)' | '(b7-7)' | '(b7-8)' | '(b7-9)'
	| '(b8-1)' | '(b8-2)' | '(b8-3)' | '(b8-4)' | '(b8-5)' | '(b8-6)' | '(b8-7)' | '(b8-8)' | '(b8-9)'
	| '(b9-1)' | '(b9-2)' | '(b9-3)' | '(b9-4)' | '(b9-5)' | '(b9-6)' | '(b9-7)' | '(b9-8)' | '(b9-9)';