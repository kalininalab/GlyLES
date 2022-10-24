import pytest

from glyles import Glycan


@pytest.fixture(scope="session")
def glycan():
    return Glycan("Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-3)[Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)"
                  "GlcNAc(b1-6)]Gal(b1-3)[GlcNAc(a1-4)Gal(b1-4)GlcNAc6S(b1-6)]GalNAc")


def test_count_match_nodes(glycan):
    assert glycan.count("Gal", match_nodes=True) == 7
    assert glycan.count(Glycan("Gal"), match_nodes=True) == 7

    assert glycan.count("GalNAc", match_nodes=True) == 7
    assert glycan.count(Glycan("GalNAc"), match_nodes=True) == 7

    assert glycan.count("Glc", match_nodes=True) == 4
    assert glycan.count(Glycan("Glc"), match_nodes=True) == 4

    assert glycan.count("GlcNAc", match_nodes=True) == 4
    assert glycan.count(Glycan("GlcNAc"), match_nodes=True) == 4

    assert glycan.count("GlcNAc6S", match_nodes=True) == 4
    assert glycan.count(Glycan("GlcNAc6S"), match_nodes=True) == 4

    assert glycan.count("Fuc", match_nodes=True) == 2
    assert glycan.count(Glycan("Fuc"), match_nodes=True) == 2


def test_count_match_root(glycan):
    assert glycan.count("Gal", match_root=True) == 1
    assert glycan.count("GalNAc", match_root=True) == 1
    assert glycan.count("Glc", match_root=True) == 0
    assert glycan.count("GlcNAc", match_root=True) == 0
    assert glycan.count("GlcNAc6S", match_root=True) == 0
    assert glycan.count("Fuc", match_root=True) == 0


def test_count_match_leaves(glycan):
    assert glycan.count("Gal", match_leaves=True) == 2
    assert glycan.count("GalNAc", match_leaves=True) == 2
    assert glycan.count("Glc", match_leaves=True) == 1
    assert glycan.count("GlcNAc", match_leaves=True) == 1
    assert glycan.count("GlcNAc6S", match_leaves=True) == 1
    assert glycan.count("Fuc", match_leaves=True) == 2


def test_count_match_nodes_all_fg(glycan):
    assert glycan.count("Gal", match_all_fg=True, match_nodes=True) == 4
    assert glycan.count("GalNAc", match_all_fg=True, match_nodes=True) == 3
    assert glycan.count("Glc", match_all_fg=True, match_nodes=True) == 0
    assert glycan.count("GlcNAc", match_all_fg=True, match_nodes=True) == 3
    assert glycan.count("GlcNAc6S", match_all_fg=True, match_nodes=True) == 1
    assert glycan.count("Fuc", match_all_fg=True, match_nodes=True) == 2


def test_count_match_root_all_fg(glycan):
    assert glycan.count("Gal", match_all_fg=True, match_root=True) == 0
    assert glycan.count("GalNAc", match_all_fg=True, match_root=True) == 1
    assert glycan.count("Glc", match_all_fg=True, match_root=True) == 0
    assert glycan.count("GlcNAc", match_all_fg=True, match_root=True) == 0
    assert glycan.count("GlcNAc6S", match_all_fg=True, match_root=True) == 0
    assert glycan.count("Fuc", match_all_fg=True, match_root=True) == 0


def test_count_match_leaves_all_fg(glycan):
    assert glycan.count("Gal", match_all_fg=True, match_leaves=True) == 0
    assert glycan.count("GalNAc", match_all_fg=True, match_leaves=True) == 2
    assert glycan.count("Glc", match_all_fg=True, match_leaves=True) == 0
    assert glycan.count("GlcNAc", match_all_fg=True, match_leaves=True) == 1
    assert glycan.count("GlcNAc6S", match_all_fg=True, match_leaves=True) == 0
    assert glycan.count("Fuc", match_all_fg=True, match_leaves=True) == 2


def test_count_match_nodes_some_fg(glycan):
    assert glycan.count("Gal", match_some_fg=True, match_nodes=True) == 7
    assert glycan.count("GalNAc", match_some_fg=True, match_nodes=True) == 3
    assert glycan.count("Glc", match_some_fg=True, match_nodes=True) == 4
    assert glycan.count("GlcNAc", match_some_fg=True, match_nodes=True) == 4
    assert glycan.count("GlcNAc6S", match_some_fg=True, match_nodes=True) == 1
    assert glycan.count("Fuc", match_some_fg=True, match_nodes=True) == 2


def test_count_match_root_some_fg(glycan):
    assert glycan.count("Gal", match_some_fg=True, match_root=True) == 1
    assert glycan.count("GalNAc", match_some_fg=True, match_root=True) == 1
    assert glycan.count("Glc", match_some_fg=True, match_root=True) == 0
    assert glycan.count("GlcNAc", match_some_fg=True, match_root=True) == 0
    assert glycan.count("GlcNAc6S", match_some_fg=True, match_root=True) == 0
    assert glycan.count("Fuc", match_some_fg=True, match_root=True) == 0


def test_count_match_leaves_some_fg(glycan):
    assert glycan.count("Gal", match_some_fg=True, match_leaves=True) == 2
    assert glycan.count("GalNAc", match_some_fg=True, match_leaves=True) == 2
    assert glycan.count("Glc", match_some_fg=True, match_leaves=True) == 1
    assert glycan.count("GlcNAc", match_some_fg=True, match_leaves=True) == 1
    assert glycan.count("GlcNAc6S", match_some_fg=True, match_leaves=True) == 0
    assert glycan.count("Fuc", match_some_fg=True, match_leaves=True) == 2


def test_count_match_nodes_poly(glycan):
    assert glycan.count("Gal(b1-4)GlcNAc", match_nodes=True) == 3
    assert glycan.count("Gal(a1-2)GlcNAc", match_nodes=True) == 3
    assert glycan.count("Gal(b1-4)GlcNAc6S", match_nodes=True) == 3
    assert glycan.count("Gal(a1-2)GlcNAc6S", match_nodes=True) == 3
    assert glycan.count("GlcNAc(b1-6)Gal", match_nodes=True) == 4
    assert glycan.count("GlcNAc(a1-2)Gal", match_nodes=True) == 4
    assert glycan.count("GlcNAc(a1-4)Gal", match_nodes=True) == 4


def test_count_match_edges(glycan):
    assert glycan.count("Gal(b1-4)GlcNAc", match_edges=True, match_nodes=True) == 3
    assert glycan.count("Gal(a1-2)GlcNAc", match_edges=True, match_nodes=True) == 0
    assert glycan.count("Gal(b1-4)GlcNAc6S", match_edges=True, match_nodes=True) == 3
    assert glycan.count("Gal(a1-2)GlcNAc6S", match_edges=True, match_nodes=True) == 0
    assert glycan.count("GlcNAc(b1-6)Gal", match_edges=True, match_nodes=True) == 2
    assert glycan.count("GlcNAc(b1-3)Gal", match_edges=True, match_nodes=True) == 1
    assert glycan.count("GlcNAc(a1-4)Gal", match_edges=True, match_nodes=True) == 1


def test_count_match_edges_all_fg(glycan):
    assert glycan.count("Gal(b1-4)GlcNAc", match_edges=True, match_all_fg=True, match_nodes=True) == 2
    assert glycan.count("Gal(a1-2)GlcNAc", match_edges=True, match_all_fg=True, match_nodes=True) == 0
    assert glycan.count("Gal(b1-4)GlcNAc6S", match_edges=True, match_all_fg=True, match_nodes=True) == 1
    assert glycan.count("Gal(a1-2)GlcNAc6S", match_edges=True, match_all_fg=True, match_nodes=True) == 0
    assert glycan.count("GlcNAc(b1-6)Gal", match_edges=True, match_all_fg=True, match_nodes=True) == 1
    assert glycan.count("GlcNAc(b1-3)Gal", match_edges=True, match_all_fg=True, match_nodes=True) == 1
    assert glycan.count("GlcNAc(a1-4)Gal", match_edges=True, match_all_fg=True, match_nodes=True) == 1


def test_count_match_edges_some_fg(glycan):
    assert glycan.count("Gal(b1-4)GlcNAc", match_edges=True, match_some_fg=True, match_nodes=True) == 3
    assert glycan.count("Gal(a1-2)GlcNAc", match_edges=True, match_some_fg=True, match_nodes=True) == 0
    assert glycan.count("Gal(b1-4)GlcNAc6S", match_edges=True, match_some_fg=True, match_nodes=True) == 1
    assert glycan.count("Gal(a1-2)GlcNAc6S", match_edges=True, match_some_fg=True, match_nodes=True) == 0
    assert glycan.count("GlcNAc(b1-6)Gal", match_edges=True, match_some_fg=True, match_nodes=True) == 2
    assert glycan.count("GlcNAc(b1-3)Gal", match_edges=True, match_some_fg=True, match_nodes=True) == 1
    assert glycan.count("GlcNAc(a1-4)Gal", match_edges=True, match_some_fg=True, match_nodes=True) == 1


def test_count_functional_groups(glycan):
    assert glycan.count_functional_groups("S") == 1
    assert glycan.count_functional_groups("Ac") == 7
    assert glycan.count_functional_groups("NC(=O)C") == 7
    assert glycan.count_functional_groups("[#6]-[#8]-[#1]") == 31


def test_count_protonation(glycan):
    assert glycan.count_protonation(groups=False) == 2
    assert glycan.count_protonation(groups=True) == 1
