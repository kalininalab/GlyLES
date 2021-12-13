from glyles.glycans.utils import Config
from glyles.grammar.parse import Glycan


def check_initial(g, name, num_children, config=None):
    assert g.nodes[0]["type"].get_name() == name
    assert len(g.edges(0)) == num_children

    if config is not None:
        assert g.nodes[0]["type"].get_config() == config


def check_child(g, id_parent, id_child, name, edge, num_children):
    assert g.get_edge_data(id_parent, id_child)["type"] == edge
    assert g.nodes[id_child]["type"].get_name() == name
    assert len(g.edges(id_child)) == num_children


def split_children(g, id_children, child_1):
    if g.nodes[id_children[0]]["type"].get_name() == child_1:
        id_child_1, id_child_2 = id_children
    else:
        id_child_2, id_child_1 = id_children
    return id_child_1, id_child_2


class TestParser:
    def test_parse_1(self):
        g = Glycan("Man").get_tree()

        check_initial(g, "Man", 0, Config.UNDEF)

    def test_parse_2(self):
        g = Glycan("Man(a1-4)Glc").get_tree()

        check_initial(g, "Glc", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Man", "(a1-4)", 0)

    def test_parse_3(self):
        g = Glycan("Man(a1-4)Glc(a1-3)Tal").get_tree()

        check_initial(g, "Tal", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Glc", "(a1-3)", 1)
        id_child_2 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_2, "Man", "(a1-4)", 0)

    def test_parse_4(self):
        g = Glycan("Man(a1-4)[Glc(a1-3)]Tal").get_tree()

        check_initial(g, "Tal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Glc")

        check_child(g, 0, id_child_1, "Glc", "(a1-3)", 0)
        check_child(g, 0, id_child_2, "Man", "(a1-4)", 0)

    def test_parse_5(self):
        g = Glycan("Man(a1-2)[Glc(a1-3)Tal(b1-4)]Gal").get_tree()

        check_initial(g, "Gal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Tal")

        check_child(g, 0, id_child_1, "Tal", "(b1-4)", 1)
        check_child(g, 0, id_child_2, "Man", "(a1-2)", 0)

        id_child_11 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_11, "Glc", "(a1-3)", 0)

    def test_parse_6(self):
        g = Glycan("Man(a1-2)Glc(a1-3)[Tal(b1-4)]Gal").get_tree()

        check_initial(g, "Gal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Tal")

        check_child(g, 0, id_child_1, "Tal", "(b1-4)", 0)
        check_child(g, 0, id_child_2, "Glc", "(a1-3)", 1)

        id_child_21 = list(g.edges(id_child_2))[0][1]
        check_child(g, id_child_2, id_child_21, "Man", "(a1-2)", 0)

    def test_parse_7(self):
        g = Glycan("Man(a1-4)Glc(a1-2)[Tal(b1-4)Gal(b1-3)]Tal").get_tree()

        check_initial(g, "Tal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Gal")

        check_child(g, 0, id_child_1, "Gal", "(b1-3)", 1)
        check_child(g, 0, id_child_2, "Glc", "(a1-2)", 1)

        id_child_11 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_11, "Tal", "(b1-4)", 0)

        id_child_21 = list(g.edges(id_child_2))[0][1]
        check_child(g, id_child_2, id_child_21, "Man", "(a1-4)", 0)

    def test_parse_8(self):
        g = Glycan("Man(a1-2)[Glc(a1-3)Tal(b1-4)]Gal(b1-3)Tal").get_tree()

        check_initial(g, "Tal", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Gal", "(b1-3)", 2)

        id_children_1 = [x[1] for x in list(g.edges(id_child_1))]
        id_child_11, id_child_12 = split_children(g, id_children_1, "Tal")

        check_child(g, id_child_1, id_child_11, "Tal", "(b1-4)", 1)
        check_child(g, id_child_1, id_child_12, "Man", "(a1-2)", 0)

        id_child_111 = list(g.edges(id_child_11))[0][1]
        check_child(g, id_child_11, id_child_111, "Glc", "(a1-3)", 0)

    def test_parse_9(self):
        g = Glycan("Man a").get_tree()

        check_initial(g, "Man", 0, Config.ALPHA)

    def test_parse_10(self):
        g = Glycan("Man b").get_tree()

        check_initial(g, "Man", 0, Config.BETA)

    def test_parse_11(self):
        g = Glycan("Man(a1-4)Glc a").get_tree()

        check_initial(g, "Glc", 1, Config.ALPHA)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Man", "(a1-4)", 0)

    def test_parse_12(self):
        g = Glycan("Man(a1-4)Glc b").get_tree()

        check_initial(g, "Glc", 1, Config.BETA)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Man", "(a1-4)", 0)

    def test_parse_13(self):
        g = Glycan("Fuc(a1-2)Gal(b1-3)GalNAc").get_tree()

        check_initial(g, "GalNAc", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Gal", "(b1-3)", 1)
        id_child_2 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_2, "Fuc", "(a1-2)", 0)

    def test_parse_14(self):
        g = Glycan("Fuc(a1-2)Gal(b1-3)GlcNAc6S").get_tree()

        check_initial(g, "GlcNAc6S", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Gal", "(b1-3)", 1)
        id_child_2 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_2, "Fuc", "(a1-2)", 0)

    def test_parse_15(self):
        g = Glycan("Fuc(a1-2)Gal(b1-4)Gal6S").get_tree()

        check_initial(g, "Gal6S", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Gal", "(b1-4)", 1)
        id_child_2 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_2, "Fuc", "(a1-2)", 0)

    def test_parse_16(self):
        g = Glycan("Fuc(a1-2)Gal(a1-3)[Fuc(a1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc").get_tree()

        check_initial(g, "GlcNAc", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "GlcNAc", "(b1-4)", 1)
        id_child_2 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_2, "Man", "(b1-4)", 2)

        id_children_2 = [x[1] for x in list(g.edges(id_child_2))]
        id_child_31, id_child_32 = split_children(g, id_children_2, "Man")

        check_child(g, id_child_2, id_child_31, "Man", "(a1-6)", 1)
        check_child(g, id_child_2, id_child_32, "Gal", "(a1-3)", 1)

        id_child_311 = list(g.edges(id_child_31))[0][1]
        check_child(g, id_child_31, id_child_311, "Fuc", "(a1-2)", 0)

        id_child_321 = list(g.edges(id_child_32))[0][1]
        check_child(g, id_child_32, id_child_321, "Fuc", "(a1-2)", 0)

    def test_parse_17(self):
        g = Glycan("Fuc(a1-2)Gal(a1-3)[Fuc(a1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc").get_tree()

        check_initial(g, "GlcNAc", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "GlcNAc")

        check_child(g, 0, id_child_1, "GlcNAc", "(b1-4)", 1)
        check_child(g, 0, id_child_2, "Fuc", "(a1-6)", 0)

        id_child_2 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_2, "Man", "(b1-4)", 2)
