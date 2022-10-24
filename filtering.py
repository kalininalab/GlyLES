import random

with open("tests/data/misc/pubchem_check_1.tsv", "r") as pubchem, open("tests/data/misc/hard_2.tsv", "w") as output:
    for line in pubchem.readlines():
        try:
            iupac, smiles, pid = line.strip().split("\t")[:3]
            if "-onic" in iupac or "-aric" in iupac or "-ulosonic" in iupac or "-ulosaric" in iupac:
                print(iupac, smiles, pid, sep="\t", file=output)
        except:
            pass
exit(0)

count = [0, 0]
with open("tests/data/general.tsv", "r") as gen, open("tests/data/pubchem_mono.tsv", "r") as mon, \
        open("tests/data/pubchem_poly.tsv", "r") as pol, open("tests/data/glycam.tsv", "r") as cam, \
        open("tests/data/anhydro.tsv", "r") as anh, open("runtime.tsv", "w") as output:
    glycans = []
    for file in [gen, mon, pol, cam, anh]:
        for line in file.readlines():
            count[0] += 1
            if random.random() < 0.2:
                count[1] += 1
                try:
                    iupac, smiles = line.strip().split("\t")[:2]
                    print(iupac, smiles, sep="\t", file=output)
                except:
                    pass

print(count)