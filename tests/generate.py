import json
from typing import Union

from glyles.grammar.parse import Glycan
from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.utils import ParseError


def filter_done(json_query, test_stream, look_stream):
    data = json.load(open(json_query, "r"))
    print("JSON read in completed")
    f = MonomerFactory()
    for i, d in enumerate(data):
        print(f"\r{i} / {len(data)}", end="")
        try:
            names = [d["cmpdname"]]
            if "cmpdsynonym" in d:
                names += d["cmpdsynonym"]
            found = False
            for n in names:
                try:
                    if Glycan(n, f, tree_only=True).get_tree() is not None:
                        print("\t".join([n, d["isosmiles"], d['cid']]), file=test_stream)
                        test_stream.flush()
                        found = True
                except ParseError:
                    pass
            if not found:
                print("\t".join([d["cmpdname"], d["isosmiles"], d['cid']]), file=look_stream)
                look_stream.flush()
        except Exception as e:
            print(f"\n{e}")
            pass


def main(mode):
    if mode == 1:
        with open("data/pubchem_mono.tsv", "r") as read:
            monomers = set([line.split("\t")[1] for line in read.readlines()])
        with open("data/pubchem_check.tsv", "r") as look:
            monomers = monomers.union(set([line.split("\t")[1] for line in look.readlines()]))
        with open("data/pubchem_poly_results.tsv", "r") as res, open("query.txt", "w") as out:
            data = {}
            for line in res.readlines():
                if "." in line:
                    continue
                parts = line.strip().split("\t")
                if parts[1] not in monomers:
                    data[parts[1]] = parts[0]
            for v in data.values():
                print(f"{v}", file=out)
    elif mode == 2:
        with open("data/pubchem_poly.tsv", "w") as read, open("data/pubchem_check.tsv", "a") as look:
            filter_done("data/PubChem_Compound_polymer_2.json", read, look)
    else:
        with open("data/PubChem_compound_polymer.json", "r") as poly, open("data/PubChem_compound_polymer_2.json", "w") as out:
            for line in poly.readlines():
                if line == "\t}\n":
                    out.write("\t},\n")
                else:
                    out.write(line)


if __name__ == '__main__':
    # main(2)
    with open("data/profiling.tsv", "w") as out:
        for filename in [
            "data/glycowork_mono.txt", "data/glycowork_poly.txt", "data/pubchem_mono.tsv", "data/pubchem_poly.tsv"
        ]:
            with open(filename, "r") as data:
                lines = sorted(data.readlines(), key=lambda l: len(l), reverse=True)
                print("".join(lines[:len(lines) // 10]), file=out)
