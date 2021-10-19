from glyles.grammar.parse import parse

smiles = [
    "Gal(a1-4)[Glc(b1-3)]Man",
    "Man(a1-4)Man",
    "Man(a1-4)Man(a1-4)Man",
    "Gal(a1-4)[Glc(b1-3)Fru(a1-4)]Man",
    "Man(a1-4)Gal(a1-4)[Glc(b1-3)]Man",
    "Man(a1-4)Gal(a1-4)[Glc(b1-3)Man(a1-4)]Man",
    "Man(a1-4)[Glc(b1-3)Man(a1-4)]Gal(a1-4)Man",
]

# parse(smiles[-1])
# exit(0)
for smile in smiles:
    print(smile)
    g = parse(smile)
