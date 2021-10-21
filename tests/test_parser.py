from glyles.grammar.parse import parse

smiles = [
    "Man",
    "Man(a1-4)Glc",
    "Man(a1-4)Glc(a1-3)Fru",
    "Man(a1-4)[Glc(a1-3)]Fru",
    "Man(a1-4)[Glc(a1-3)Fru(b1-4)]Gal",
    "Man(a1-4)Glc(a1-3)[Fru(b1-4)]Gal",
    "Man(a1-4)Glc(a1-3)[Fru(b1-4)Gal(b1-3)]Xyl",
    "Man(a1-4)[Glc(a1-3)Fru(b1-4)]Gal(b1-3)Xyl",
]


def check_initial(g, name, num_children):
    assert g.nodes[0]["type"].value["name"] == name
    assert len(g.edges(0)) == num_children


def check_child(g, id_parent, id_child, name, edge, num_children):
    assert g.get_edge_data(id_parent, id_child)["type"] == edge
    assert g.nodes[id_child]["type"].value["name"] == name
    assert len(g.edges(id_child)) == num_children


def split_children(g, id_children, child_1):
    if g.nodes[id_children[0]]["type"].value["name"] == child_1:
        id_child_1, id_child_2 = id_children
    else:
        id_child_2, id_child_1 = id_children
    return id_child_1, id_child_2


class TestParser:
    def test_parse_1(self):
        g = parse(smiles[0])

        check_initial(g, "Man", 0)

    def test_parse_2(self):
        g = parse(smiles[1])

        check_initial(g, "Glc", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Man", "(a1-4)", 0)

    def test_parse_3(self):
        g = parse(smiles[2])

        check_initial(g, "Fru", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Glc", "(a1-3)", 1)
        id_child_2 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_2, "Man", "(a1-4)", 0)

    def test_parse_4(self):
        g = parse(smiles[3])

        check_initial(g, "Fru", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Glc")

        check_child(g, 0, id_child_1, "Glc", "(a1-3)", 0)
        check_child(g, 0, id_child_2, "Man", "(a1-4)", 0)

    def test_parse_5(self):
        g = parse(smiles[4])

        check_initial(g, "Gal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Fru")

        check_child(g, 0, id_child_1, "Fru", "(b1-4)", 1)
        check_child(g, 0, id_child_2, "Man", "(a1-4)", 0)

        id_child_11 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_11, "Glc", "(a1-3)", 0)

    def test_parse_6(self):
        g = parse(smiles[5])

        check_initial(g, "Gal", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Fru")

        check_child(g, 0, id_child_1, "Fru", "(b1-4)", 0)
        check_child(g, 0, id_child_2, "Glc", "(a1-3)", 1)

        id_child_21 = list(g.edges(id_child_2))[0][1]
        check_child(g, id_child_2, id_child_21, "Man", "(a1-4)", 0)

    def test_parse_7(self):
        g = parse(smiles[6])

        check_initial(g, "Xyl", 2)
        id_children_1 = [x[1] for x in list(g.edges(0))]
        id_child_1, id_child_2 = split_children(g, id_children_1, "Gal")

        check_child(g, 0, id_child_1, "Gal", "(b1-3)", 1)
        check_child(g, 0, id_child_2, "Glc", "(a1-3)", 1)

        id_child_11 = list(g.edges(id_child_1))[0][1]
        check_child(g, id_child_1, id_child_11, "Fru", "(b1-4)", 0)

        id_child_21 = list(g.edges(id_child_2))[0][1]
        check_child(g, id_child_2, id_child_21, "Man", "(a1-4)", 0)

    def test_parse_8(self):
        g = parse(smiles[7])

        check_initial(g, "Xyl", 1)
        id_child_1 = list(g.edges(0))[0][1]
        check_child(g, 0, id_child_1, "Gal", "(b1-3)", 2)

        id_children_1 = [x[1] for x in list(g.edges(id_child_1))]
        id_child_11, id_child_12 = split_children(g, id_children_1, "Fru")

        check_child(g, id_child_1, id_child_11, "Fru", "(b1-4)", 1)
        check_child(g, id_child_1, id_child_12, "Man", "(a1-4)", 0)

        id_child_111 = list(g.edges(id_child_11))[0][1]
        check_child(g, id_child_11, id_child_111, "Glc", "(a1-3)", 0)
