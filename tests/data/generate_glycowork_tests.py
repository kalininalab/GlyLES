import os
import pickle
import re


def check_parenthesis(expression):
    open_tup = tuple('({[')
    close_tup = tuple(')}]')
    m = dict(zip(open_tup, close_tup))
    queue = []

    for i in expression:
        if i in open_tup:
            queue.append(m[i])
        elif i in close_tup:
            if not queue or i != queue.pop():
                return False
    if not queue:
        return True
    return False


def check_iupac(iupac):
    return iupac[0] != "[" and check_parenthesis(iupac) and iupac != "dHex" \
        and all(x not in iupac for x in [
            "OMe", "OAc", "OP", "bond", "0dHex", "(z", "-z", "GalAP3", "Galc",
            "5)Ara", "5)]Ara", "5)D-Ara", "8)Ara", "6)Fru", ")5Rib"
            "DL-", "Gal6ulo",

            "1)Glc(a1",
            "2)[FucNAc(a1-2)]", "2)[Sug(b1-3)Gal(a1-2)]", "2)[Glc(b1-2)]", "2)[Xyl(b1-2)]", "2)[Man(a1-2)Man(a1-2)]",
            "3)[Col(b1-3)]", "3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)]", "3)[D-4dAraHex(b1-3)]", "3)[Neu5Ac(a2-3)]",
            "3)[Gal(b1-3)GlcNAc(b1-3)]", "3)[Gal(b1-4)GlcNAc(b1-3)]", "[Rha(a1-3)][Gal(b1-3)]", "3)[Neu5Gc(a2-3)]",
            "3)[Man(a1-2)Man(a1-3)]", "3)[Glc(b1-3)]", "3)[Gal(b1-3)]", "3)[Man(a1-3)]",
            "4)[Fuc(a1-4)]", "4)[GalA(a1-4)]", "4)[Rha(b1-4)]", "4)[Gal(a1-4)]",
            "[GlcPGro(b1-6)][Gal(b1-6)]", "6)[Glc(b1-6)]", "6)[Man(a1-6)Man(a1-6)]", "6)[Man(a1-2)Man(a1-6)]",
            "6)[Man(a1-2)Man(a1-2)Man(a1-2)Man(a1-6)]", "6)[Xyl(a1-6)]",
        ]) \
        and not iupac[-1].isdigit() and len(re.findall(r"Glc[a-zA-Z\d]*1[a-zA-Z]*[a-zA-Z\d]*", iupac)) == 0 \
        and len(re.findall(r"Qui[a-zA-Z\d]*Ac?", iupac)) == 0 \
        and len(re.findall(r"Ara[a-zA-Z\d]*Ac?", iupac)) == 0 \
        and len(re.findall(r"Fuc[a-zA-Z\d]*Ac?", iupac)) == 0 \
        and len(re.findall(r"Rha[a-zA-Z\d]*Ac?", iupac)) == 0 \
        and len(re.findall(r"1,\d-Anhydro-", iupac)) == 0


def load_monomers(filepath):
    os.system("wget https://github.com/BojarLab/glycowork/raw/master/glycowork/glycan_data/lib_v6.pkl")
    with open(filepath, "w") as output:
        for iupac in pickle.load(open("lib_v6.pkl", "rb")):
            if check_iupac(iupac):
                output.write(f"{iupac}\n")
    os.system("rm lib_v6.pkl")


def load_polymers(filepath):
    os.system("wget https://raw.githubusercontent.com/BojarLab/glycowork/master/glycowork/glycan_data/v6_sugarbase.csv")
    with open("v6_sugarbase.csv", "r") as data, open(filepath, "w") as output:
        for line in data.readlines()[1:]:
            iupac = re.findall("\".*?\"", line)
            if len(iupac) > 0:
                iupac = iupac[0][1:-1]
            else:
                iupac = line.split(",")[0]
            if check_iupac(iupac):
                output.write(f"{iupac}\n")
    os.system("rm v6_sugarbase.csv")


load_monomers("misc/glycowork_mono_test.txt")
load_polymers("misc/glycowork_poly_test.txt")
