import numpy as np
import pytest

from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.utils import Config, Lactole
from glyles import Glycan
from tests.utils import reduce_notation


def check_initial(g, name, num_children, config=None, lactole=None):
    assert g.nodes[0]["type"].get_name() == name
    assert len(g.edges(0)) == num_children

    if config is not None:
        assert g.nodes[0]["type"].get_config() == config
    if lactole is not None and name != "Api":
        assert g.nodes[0]["type"].get_lactole() == lactole


def check_child(g, id_parent, id_child, name, edge, num_children, lactole=None):
    assert g.get_edge_data(id_parent, id_child)["type"] == edge
    assert g.nodes[id_child]["type"].get_name() == name
    assert len(g.edges(id_child)) == num_children

    if lactole is not None and name != "Api":
        assert g.nodes[id_child]["type"].get_lactole() == lactole


def split_children(g, id_children, child_1):
    if g.nodes[id_children[0]]["type"].get_name() == child_1:
        id_child_1, id_child_2 = id_children
    else:
        id_child_2, id_child_1 = id_children
    return id_child_1, id_child_2


def split_ternary_children(g, id_children, child_1, child_2):
    if g.nodes[id_children[0]]["type"].get_name() == child_1:
        if g.nodes[id_children[1]]["type"].get_name() == child_2:
            id_child_1, id_child_2, id_child_3 = id_children
        else:
            id_child_1, id_child_3, id_child_2 = id_children
    elif g.nodes[id_children[1]]["type"].get_name() == child_1:
        if g.nodes[id_children[0]]["type"].get_name() == child_2:
            id_child_2, id_child_1, id_child_3 = id_children
        else:
            id_child_3, id_child_1, id_child_2 = id_children
    else:
        if g.nodes[id_children[0]]["type"].get_name() == child_2:
            id_child_2, id_child_3, id_child_1 = id_children
        else:
            id_child_3, id_child_2, id_child_1 = id_children
    return id_child_1, id_child_2, id_child_3


def test_grammar_simple():
    g = Glycan("3dMan5S6Pfa").get_tree()

    check_initial(g, "Man", 0, Config.ALPHA, lactole=Lactole.FURANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_extreme_1(mode):
    g = Glycan(reduce_notation("Gal(a1-4)" * 9 + "Glc", mode))
    check_initial(g.get_tree(), "Glc", 1, Config.UNDEF, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_extreme_2(mode):
    g = Glycan(reduce_notation("Gal(a1-4)" * 8 + "[" + "Man(b1-3)" * 8 + "]" + "Glc", mode))
    check_initial(g.get_tree(), "Glc", 2, Config.UNDEF, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("iupac", [
    "Man", "Man a", "Man b", "Mana", "Manb", "Manp", "Manf", "Manpa", "Manpb", "Manfa", "Manfb", "Man(a1-2)Gal",
    "Man(a1-2)Gal a", "Man(a1-2)Gal b", "Mana(a1-2)Gal", "Manb(a1-2)Gal", "Manp(a1-2)Gal", "Manf(a1-2)Gal",
    "Manpa(a1-2)Gal", "Manpb(a1-2)Gal", "Manfa(a1-2)Gal", "Manfb(a1-2)Gal", "Man(b1-2)Gala", "Man(b1-2)Galb",
    "Man(b1-2)Galp", "Man(b1-2)Galf", "Man(b1-2)Galpa", "Man(b1-2)Galpb", "Man(b1-2)Galfa", "Man(b1-2)Galfb"
])
@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_grammar_main(iupac, mode):
    assert Glycan(reduce_notation(iupac, mode), tree_only=True).get_tree() is not None

def test_parse_1():
    g = Glycan("Man").get_tree()

    check_initial(g, "Man", 0, Config.UNDEF, lactole=Lactole.PYRANOSE)

def test_parse_1_2():
    g = Glycan("Manpa").get_tree()

    check_initial(g, "Man", 0, Config.ALPHA, lactole=Lactole.PYRANOSE)

@pytest.mark.slow
@pytest.mark.parametrize("mono", list(MonomerFactory().pyranoses()))
@pytest.mark.parametrize("config", [Config.ALPHA, Config.BETA, Config.UNDEF])
@pytest.mark.parametrize("lactole", ["", "p"])
def test_parse_1_multi(mono, config, lactole):
    iupac = mono + lactole

    if config == Config.ALPHA:
        iupac += " a"
    elif config == Config.BETA:
        iupac += " b"
    g = Glycan(iupac, tree_only=True).get_tree()

    check_initial(g, mono, 0, config, lactole=Lactole.PYRANOSE)

@pytest.mark.slow
@pytest.mark.parametrize("mono", list(MonomerFactory().furanoses()))
@pytest.mark.parametrize("config", [Config.ALPHA, Config.BETA, Config.UNDEF])
def test_parse_1_furanose(mono, config):
    iupac = mono + "f"

    if config == Config.ALPHA:
        iupac += " a"
    elif config == Config.BETA:
        iupac += " b"
    g = Glycan(iupac, tree_only=True).get_tree()

    check_initial(g, mono, 0, config, lactole=Lactole.FURANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_2(mode):
    iupac = "Man(a1-4)Glc"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Glc", 1, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Man", "(a1-4)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_3(mode):
    iupac = "Man(a1-4)Glc(a1-3)Tal"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Tal", 1, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Glc", "(a1-3)", 1, lactole=Lactole.PYRANOSE)
    id_child_2 = list(g.edges(id_child_1))[0][1]
    check_child(g, id_child_1, id_child_2, "Man", "(a1-4)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_4(mode):
    iupac = "Man(a1-4)[Glc(a1-3)]Tal"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Tal", 2, lactole=Lactole.PYRANOSE)
    id_children_1 = [x[1] for x in list(g.edges(0))]
    id_child_1, id_child_2 = split_children(g, id_children_1, "Glc")

    check_child(g, 0, id_child_1, "Glc", "(a1-3)", 0, lactole=Lactole.PYRANOSE)
    check_child(g, 0, id_child_2, "Man", "(a1-4)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_5(mode):
    iupac = "Man(a1-2)[Glc(a1-3)Tal(b1-4)]Gal"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Gal", 2, lactole=Lactole.PYRANOSE)
    id_children_1 = [x[1] for x in list(g.edges(0))]
    id_child_1, id_child_2 = split_children(g, id_children_1, "Tal")

    check_child(g, 0, id_child_1, "Tal", "(b1-4)", 1, lactole=Lactole.PYRANOSE)
    check_child(g, 0, id_child_2, "Man", "(a1-2)", 0, lactole=Lactole.PYRANOSE)

    id_child_11 = list(g.edges(id_child_1))[0][1]
    check_child(g, id_child_1, id_child_11, "Glc", "(a1-3)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_6(mode):
    iupac = "Man(a1-2)Glc(a1-3)[Tal(b1-4)]Gal"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Gal", 2, lactole=Lactole.PYRANOSE)
    id_children_1 = [x[1] for x in list(g.edges(0))]
    id_child_1, id_child_2 = split_children(g, id_children_1, "Tal")

    check_child(g, 0, id_child_1, "Tal", "(b1-4)", 0, lactole=Lactole.PYRANOSE)
    check_child(g, 0, id_child_2, "Glc", "(a1-3)", 1, lactole=Lactole.PYRANOSE)

    id_child_21 = list(g.edges(id_child_2))[0][1]
    check_child(g, id_child_2, id_child_21, "Man", "(a1-2)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_7(mode):
    iupac = "Man(a1-4)Glc(a1-2)[Tal(b1-4)Gal(b1-3)]Tal"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Tal", 2, lactole=Lactole.PYRANOSE)
    id_children_1 = [x[1] for x in list(g.edges(0))]
    id_child_1, id_child_2 = split_children(g, id_children_1, "Gal")

    check_child(g, 0, id_child_1, "Gal", "(b1-3)", 1, lactole=Lactole.PYRANOSE)
    check_child(g, 0, id_child_2, "Glc", "(a1-2)", 1, lactole=Lactole.PYRANOSE)

    id_child_11 = list(g.edges(id_child_1))[0][1]
    check_child(g, id_child_1, id_child_11, "Tal", "(b1-4)", 0, lactole=Lactole.PYRANOSE)

    id_child_21 = list(g.edges(id_child_2))[0][1]
    check_child(g, id_child_2, id_child_21, "Man", "(a1-4)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_8(mode):
    iupac = "Man(a1-2)[Glc(a1-3)Tal(b1-4)]Gal(b1-3)Tal"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Tal", 1, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Gal", "(b1-3)", 2, lactole=Lactole.PYRANOSE)

    id_children_1 = [x[1] for x in list(g.edges(id_child_1))]
    id_child_11, id_child_12 = split_children(g, id_children_1, "Tal")

    check_child(g, id_child_1, id_child_11, "Tal", "(b1-4)", 1, lactole=Lactole.PYRANOSE)
    check_child(g, id_child_1, id_child_12, "Man", "(a1-2)", 0, lactole=Lactole.PYRANOSE)

    id_child_111 = list(g.edges(id_child_11))[0][1]
    check_child(g, id_child_11, id_child_111, "Glc", "(a1-3)", 0, lactole=Lactole.PYRANOSE)

def test_parse_9():
    iupac = "Man a"
    g = Glycan(iupac).get_tree()

    check_initial(g, "Man", 0, Config.ALPHA, lactole=Lactole.PYRANOSE)

def test_parse_10():
    g = Glycan("Man b").get_tree()

    check_initial(g, "Man", 0, Config.BETA, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_11(mode):
    iupac = "Man(a1-4)Glc a"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Glc", 1, Config.ALPHA, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Man", "(a1-4)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_12(mode):
    iupac = "Man(a1-4)Glc b"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Glc", 1, Config.BETA, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Man", "(a1-4)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_13(mode):
    iupac = "Fuc(a1-2)Gal(b1-3)GalNAc"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Gal", 1, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Gal", "(b1-3)", 1, lactole=Lactole.PYRANOSE)
    id_child_2 = list(g.edges(id_child_1))[0][1]
    check_child(g, id_child_1, id_child_2, "Fuc", "(a1-2)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_14(mode):
    iupac = "Fuc(a1-2)Gal(b1-3)GlcNAc6S"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Glc", 1, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Gal", "(b1-3)", 1, lactole=Lactole.PYRANOSE)
    id_child_2 = list(g.edges(id_child_1))[0][1]
    check_child(g, id_child_1, id_child_2, "Fuc", "(a1-2)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_15(mode):
    iupac = "Fuc(a1-2)Gal(b1-4)Gal6S"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Gal", 1, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Gal", "(b1-4)", 1, lactole=Lactole.PYRANOSE)
    id_child_2 = list(g.edges(id_child_1))[0][1]
    check_child(g, id_child_1, id_child_2, "Fuc", "(a1-2)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_16(mode):
    iupac = "Fuc(a1-2)Gal(a1-3)[Fuc(a1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Glc", 1, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Glc", "(b1-4)", 1, lactole=Lactole.PYRANOSE)
    id_child_2 = list(g.edges(id_child_1))[0][1]
    check_child(g, id_child_1, id_child_2, "Man", "(b1-4)", 2, lactole=Lactole.PYRANOSE)

    id_children_2 = [x[1] for x in list(g.edges(id_child_2))]
    id_child_31, id_child_32 = split_children(g, id_children_2, "Man")

    check_child(g, id_child_2, id_child_31, "Man", "(a1-6)", 1, lactole=Lactole.PYRANOSE)
    check_child(g, id_child_2, id_child_32, "Gal", "(a1-3)", 1, lactole=Lactole.PYRANOSE)

    id_child_311 = list(g.edges(id_child_31))[0][1]
    check_child(g, id_child_31, id_child_311, "Fuc", "(a1-2)", 0, lactole=Lactole.PYRANOSE)

    id_child_321 = list(g.edges(id_child_32))[0][1]
    check_child(g, id_child_32, id_child_321, "Fuc", "(a1-2)", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_17(mode):
    iupac = "Fuc(a1-2)Gal(a1-3)[Fuc(a1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Glc", 2, lactole=Lactole.PYRANOSE)
    id_children_1 = [x[1] for x in list(g.edges(0))]
    id_child_1, id_child_2 = split_children(g, id_children_1, "Glc")

    check_child(g, 0, id_child_1, "Glc", "(b1-4)", 1, lactole=Lactole.PYRANOSE)
    check_child(g, 0, id_child_2, "Fuc", "(a1-6)", 0, lactole=Lactole.PYRANOSE)

    id_child_2 = list(g.edges(id_child_1))[0][1]
    check_child(g, id_child_1, id_child_2, "Man", "(b1-4)", 2, lactole=Lactole.PYRANOSE)

@pytest.mark.fuzzy
@pytest.mark.parametrize(
    "monomers", np.random.choice(list(MonomerFactory().monomer_names()), size=500).reshape(100, 5)
)
@pytest.mark.parametrize("orientation", [Config.ALPHA, Config.BETA, Config.UNDEF])
@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_fuzzy(monomers, orientation, mode):
    c = ["a1-4", "a1-4", "a1-3", "a1-4"]

    iupac = f"{monomers[0]}({c[0]})[{monomers[1]}({c[1]}){monomers[2]}({c[2]})]{monomers[3]}({c[3]}){monomers[4]}"
    if orientation == Config.ALPHA:
        iupac += " a"
    elif orientation == Config.BETA:
        iupac += " b"

    g = Glycan(reduce_notation(iupac, mode), tree_only=True).get_tree()

    check_initial(g, monomers[4], 1, orientation, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, monomers[3], f"({c[3]})", 2, lactole=Lactole.PYRANOSE)

    id_children_1 = [x[1] for x in list(g.edges(id_child_1))]
    id_child_11, id_child_12 = split_children(g, id_children_1, monomers[2])

    check_child(g, id_child_1, id_child_11, monomers[2], f"({c[2]})", 1, lactole=Lactole.PYRANOSE)
    check_child(g, id_child_1, id_child_12, monomers[0], f"({c[0]})", 0, lactole=Lactole.PYRANOSE)

    id_child_111 = list(g.edges(id_child_11))[0][1]
    check_child(g, id_child_11, id_child_111, monomers[1], f"({c[1]})", 0, lactole=Lactole.PYRANOSE)

def test_parse_fuzzy_check():
    monomers = ['Tag', 'Mur', 'Hex', 'Kdo', 'Dig']
    orientation = Config.ALPHA
    mode = "simple"
    c = ["a1-4", "a1-4", "a1-3", "a1-4"]

    iupac = f"{monomers[0]}({c[0]})[{monomers[1]}({c[1]}){monomers[2]}({c[2]})]{monomers[3]}({c[3]}){monomers[4]}"
    if orientation == Config.ALPHA:
        iupac += " a"
    elif orientation == Config.BETA:
        iupac += " b"

    g = Glycan(reduce_notation(iupac, mode), tree_only=True).get_tree()

    check_initial(g, monomers[4], 1, orientation, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, monomers[3], f"({c[3]})", 2, lactole=Lactole.PYRANOSE)

    id_children_1 = [x[1] for x in list(g.edges(id_child_1))]
    id_child_11, id_child_12 = split_children(g, id_children_1, monomers[2])

    check_child(g, id_child_1, id_child_11, monomers[2], f"({c[2]})", 1, lactole=Lactole.PYRANOSE)
    check_child(g, id_child_1, id_child_12, monomers[0], f"({c[0]})", 0, lactole=Lactole.PYRANOSE)

    id_child_111 = list(g.edges(id_child_11))[0][1]
    check_child(g, id_child_11, id_child_111, monomers[1], f"({c[1]})", 0, lactole=Lactole.PYRANOSE)

@pytest.mark.parametrize("orientation", [Config.ALPHA, Config.BETA, Config.UNDEF])
@pytest.mark.parametrize("pos_man", [1, 2, 3, 4, 6])
@pytest.mark.parametrize("pos_glc", [2, 3, 4, 6])
@pytest.mark.parametrize("conf_glc", [Config.ALPHA, Config.BETA, Config.UNDEF])
@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_connections(orientation, pos_man, pos_glc, conf_glc, mode):
    if orientation == Config.ALPHA:
        config = "a"
    elif orientation == Config.BETA:
        config = "b"
    else:
        if mode != "full":
            return
        config = ""
    iupac = f"Man({config}{pos_man}-{pos_glc})Glc"
    if conf_glc == Config.ALPHA:
        iupac += " a"
    elif conf_glc == Config.BETA:
        iupac += " b"

    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Glc", 1, conf_glc, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Man", f"({config}{pos_man}-{pos_glc})", 0)

@pytest.mark.fuzzy
@pytest.mark.parametrize("lactoles", list(zip(*(
        np.random.choice(list(MonomerFactory().pyranose_names()), size=500).reshape(100, 5),
        np.random.choice(list(MonomerFactory().furanose_names()), size=500).reshape(100, 5),
        np.random.choice([Lactole.PYRANOSE, Lactole.FURANOSE], size=500).reshape(100, 5),
))))
@pytest.mark.parametrize("orientation", [Config.ALPHA, Config.BETA, Config.UNDEF])
@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_fuzzy_pyranoses(lactoles, orientation, mode):
    c = ["a1-4", "a1-4", "a1-3", "a1-4"]

    monomers = [lactoles[0][i] + "p" if lactoles[2][i] == Lactole.PYRANOSE else lactoles[1][i] + "f" for i in
                range(len(lactoles[0]))]

    iupac = f"{monomers[0]}({c[0]})[{monomers[1]}({c[1]}){monomers[2]}({c[2]})]{monomers[3]}({c[3]}){monomers[4]}"
    if orientation == Config.ALPHA:
        iupac += " a"
    elif orientation == Config.BETA:
        iupac += " b"

    g = Glycan(reduce_notation(iupac, mode), tree_only=True).get_tree()

    check_initial(g, monomers[4][:-1], 1, orientation, lactole=lactoles[2][4])
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, monomers[3][:-1], f"({c[3]})", 2, lactole=lactoles[2][3])

    id_children_1 = [x[1] for x in list(g.edges(id_child_1))]
    id_child_11, id_child_12 = split_children(g, id_children_1, monomers[2][:-1])

    check_child(g, id_child_1, id_child_11, monomers[2][:-1], f"({c[2]})", 1, lactole=lactoles[2][2])
    check_child(g, id_child_1, id_child_12, monomers[0][:-1], f"({c[0]})", 0, lactole=lactoles[2][0])

    id_child_111 = list(g.edges(id_child_11))[0][1]
    check_child(g, id_child_11, id_child_111, monomers[1][:-1], f"({c[1]})", 0, lactole=lactoles[2][1])

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_ternary_branching_1(mode):
    iupac = "Alt(a1-2)[Glc(a1-4)][Gal(a1-6)]Man"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Man", 3, lactole=Lactole.PYRANOSE)

    id_children_1 = [x[1] for x in list(g.edges(0))]
    id_child_11, id_child_12, id_child_13 = split_ternary_children(g, id_children_1, "Alt", "Glc")
    check_child(g, 0, id_child_11, "Alt", "(a1-2)", 0, Lactole.PYRANOSE)
    check_child(g, 0, id_child_12, "Glc", "(a1-4)", 0, Lactole.PYRANOSE)
    check_child(g, 0, id_child_13, "Gal", "(a1-6)", 0, Lactole.PYRANOSE)

@pytest.mark.parametrize("mode", ["full", "condensed", "simple"])
def test_parse_ternary_branching_2(mode):
    iupac = "Alt(a1-2)[Glc(a1-4)][Gal(a1-6)]Gul(a1-4)Man"
    g = Glycan(reduce_notation(iupac, mode)).get_tree()

    check_initial(g, "Man", 1, lactole=Lactole.PYRANOSE)
    id_child_1 = list(g.edges(0))[0][1]
    check_child(g, 0, id_child_1, "Gul", "(a1-4)", 3, lactole=Lactole.PYRANOSE)
    id_children_2 = [x[1] for x in list(g.edges(id_child_1))]
    id_child_21, id_child_22, id_child_23 = split_ternary_children(g, id_children_2, "Alt", "Glc")
    check_child(g, id_child_1, id_child_21, "Alt", "(a1-2)", 0, Lactole.PYRANOSE)
    check_child(g, id_child_1, id_child_22, "Glc", "(a1-4)", 0, Lactole.PYRANOSE)
    check_child(g, id_child_1, id_child_23, "Gal", "(a1-6)", 0, Lactole.PYRANOSE)

def test_parsing_error(caplog):
    iupac = "Alt(a1-2)[Glc(a1-4)][Gal(a1-6)]Gul(a1-4)M*#$s'\\d ;«]as;an"  # Invalid IUPAC string!
    g = Glycan(iupac).get_tree()
    assert g is None
    assert "A parsing error occurred" in caplog.records[0].msg
