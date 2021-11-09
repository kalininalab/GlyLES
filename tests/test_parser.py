from glyles.grammar.parse import Glycan

smiles = [
    "Man",
    "Man(a1-4)Glc",
    "Man(a1-4)Glc(a1-3)Tal",
    "Man(a1-4)[Glc(a1-3)]Tal",
    "Man(a1-4)[Glc(a1-3)Tal(b1-4)]Gal",
    "Man(a1-4)Glc(a1-3)[Tal(b1-4)]Gal",
    "Man(a1-4)Glc(a1-3)[Tal(b1-4)Gal(b1-3)]Tal",
    "Man(a1-4)[Glc(a1-3)Tal(b1-4)]Gal(b1-3)Tal",
]


def check_initial(g, name, num_children):
    assert g.nodes[0]["type"].get_name() == name
    assert len(g.edges(0)) == num_children


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
        g = Glycan(smiles[0], parse=False).get_tree()

        check_initial(g, "Man", 0)

    def test_parse_2(self):
        g = Glycan(smiles[1], parse=False).get_tree()

        check_initial(g, "Glc", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Man", "(a1-4)", 0)

    def test_parse_3(self):
        g = Glycan(smiles[2], parse=False).get_tree()

        check_initial(g, "Tal", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Glc", "(a1-3)", 1)
        id_child_2 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_2, "Man", "(a1-4)", 0)

    def test_parse_4(self):
        g = Glycan(smiles[3], parse=False).get_tree()

        check_initial(g, "Tal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Glc")

        check_child(g, 0, id_child_1, "Glc", "(a1-3)", 0)
        check_child(g, 0, id_child_2, "Man", "(a1-4)", 0)

    def test_parse_5(self):
        g = Glycan(smiles[4], parse=False).get_tree()

        check_initial(g, "Gal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Tal")

        check_child(g, 0, id_child_1, "Tal", "(b1-4)", 1)
        check_child(g, 0, id_child_2, "Man", "(a1-4)", 0)

        id_child_11 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_11, "Glc", "(a1-3)", 0)

    def test_parse_6(self):
        g = Glycan(smiles[5], parse=False).get_tree()

        check_initial(g, "Gal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Tal")

        check_child(g, 0, id_child_1, "Tal", "(b1-4)", 0)
        check_child(g, 0, id_child_2, "Glc", "(a1-3)", 1)

        id_child_21 = list(g.edges(id_child_2))[0][1]
        check_child(g, id_child_2, id_child_21, "Man", "(a1-4)", 0)

    def test_parse_7(self):
        g = Glycan(smiles[6], parse=False).get_tree()

        check_initial(g, "Tal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Gal")

        check_child(g, 0, id_child_1, "Gal", "(b1-3)", 1)
        check_child(g, 0, id_child_2, "Glc", "(a1-3)", 1)

        id_child_11 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_11, "Tal", "(b1-4)", 0)

        id_child_21 = list(g.edges(id_child_2))[0][1]
        check_child(g, id_child_2, id_child_21, "Man", "(a1-4)", 0)

    def test_parse_8(self):
        g = Glycan(smiles[7], parse=False).get_tree()

        check_initial(g, "Tal", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Gal", "(b1-3)", 2)

        id_children_1 = [x[1] for x in list(g.edges(id_child_1))]
        id_child_11, id_child_12 = split_children(g, id_children_1, "Tal")

        check_child(g, id_child_1, id_child_11, "Tal", "(b1-4)", 1)
        check_child(g, id_child_1, id_child_12, "Man", "(a1-4)", 0)

        id_child_111 = list(g.edges(id_child_11))[0][1]
        check_child(g, id_child_11, id_child_111, "Glc", "(a1-3)", 0)
