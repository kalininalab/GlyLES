import sys
from collections import Counter
from typing import Union, List

import networkx as nx
from PIL import ImageDraw, Image
from networkx.algorithms.isomorphism import DiGraphMatcher
import pydot
from antlr4 import InputStream, CommonTokenStream
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.mono.reactor import functional_groups
from glyles.glycans.mono.monomer import Monomer
from glyles.glycans.poly.merger import Merger
from glyles.glycans.poly.viz import Tree
from glyles.glycans.poly.walker import TreeWalker
from glyles.glycans.utils import ParseError, find_isomorphism_nx
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.grammar.GlycanParser import GlycanParser


def compare_smiles(
        c: Chem.Mol,
        s: Chem.Mol
):
    """
    Compare two molecules if they are equal.

    Args:
        c: molecule 1
        s: molecule 2

    Returns:
        True if their kekulized, canonical SMILES string is equal; False otherwise
    """
    Chem.Kekulize(c)
    Chem.Kekulize(s)

    ssmiles = Chem.MolToSmiles(s, kekuleSmiles=True)
    csmiles = Chem.MolToSmiles(c, kekuleSmiles=True)
    return csmiles == ssmiles


def recipe_equality(
        glycan: Monomer,
        query: Monomer,
        no: bool = False,
        some: bool = False,
        every: bool = False
):
    """
    Checking if two monomers are considered isomorphic under the given mode of isomorphism. One of the bool-flags has
    to be set.

    Note:
        node-match's first argument is always from first node, second argument always from second graph

    Args:
        glycan: monomer to search in
        query: monomer to find in glycan
        no: If True, only match basic monosaccharide, no functional groups
        some: If True, match all of query's groups but not all of glycan's groups
        every: If True, glycans must match exactly, i.e., modifications have to be the same

    Returns:
        True if glycan and query are considered equal under given circumstances
    """
    # check if exactly one isomorphism mode is activated
    if sum([no, some, every]) != 1:
        raise ValueError("Exactly one of arguments no, some, every has to be set to True.")

    if no:
        # if functional groups should not be considered for isomorphism, just check the root monomer of both nodes
        recipe_glycan, recipe_query = glycan.get_recipe(), query.get_recipe()
        return recipe_glycan[list(zip(*recipe_glycan))[1].index(GlycanLexer.SAC)] == \
               recipe_query[list(zip(*recipe_query))[1].index(GlycanLexer.SAC)]
    if some:
        # raise NotImplementedError("Matching to a subset of functional groups is not yet implemented. Coming soon.")
        """# if query's fgs have to be a subset of glycan's fgs, do some magic
        recipe_glycan, recipe_query = glycan.get_recipe(), query.get_recipe()
        if recipe_glycan[list(zip(*recipe_glycan))[1].index(GlycanLexer.SAC)] != \
                recipe_query[list(zip(*recipe_query))[1].index(GlycanLexer.SAC)]:
            return False
        iso = find_isomorphism_nx(glycan, query, query.name, query.c1_find)
        inv = dict([(v, k) for k, v in iso])
        ring_o = query.x[querx.x[:, 1] == 100]
        if ring_o not in inv:
            return False"""
        return len(set(query.get_recipe()).difference(set(glycan.get_recipe()))) == 0
    if every:
        # if all functional groups have to be matched, compare the smiles strings of both
        return compare_smiles(glycan.get_structure(), query.get_structure())
    return False


class Glycan:
    """
    This class is like an interaction with the Parser for the IUPAC representation of the glycan. The grammar for
    glycans is defined using ANTLR (https://www.antlr.org/). From this ANTLR is able to generate lexer and parser that
    fit the defined grammar. Don't touch those files those are auto generated and therefore mostly uncommented.

    The defined grammar discards the last glycan which is used to define the root of the glycan tree. Therefore, the
    resulting abstract syntax trees (AST)s are not intuitive.
    """

    def __init__(self, iupac, root_orientation="n", start=100, tree_only=False, full=True):
        """
        Initialize the glycan from the IUPAC string.

        Args:
            iupac (str): IUPAC string representation of the glycan to represent
            root_orientation (str): orientation of the root monomer in the glycan (choose from 'a', 'b', 'n')
            start (int): ID of the atom to start with in the root monomer when generating the SMILES
            tree_only (int): Flag indicating to only parse the tree of glycans and not the modifications
            full (bool): Flag indicating that only fully convertible glycans should be returned, i.e. all modifications
                such as 3-Anhydro-[...] are also present in the SMILES
        """
        self.iupac = iupac
        self.parse_tree = None
        self.grammar_tree = None
        self.glycan_smiles = None
        self.root_orientation = root_orientation
        self.start = start
        self.tree_only = tree_only
        self.factory = MonomerFactory()
        self.full = full
        self.tree_full = True
        self.__parse()

    def summary(self):
        """
        Aggregate some statistics of the glycan. This includes in the following order [the key in the output dictionary
        in brackets]:
        molecular formula [formula], number of atoms [atoms], number of bonds [bonds], number of rings [rings], number
        of monomers [monomers], max depth of the tree [depth], the root monomer [root], list of all leaf monomers
        [leaves], molecular weight [weight].

        Returns:
            The above named statistics are returned as dictionary with the given keys.
        """
        # generate smiles for this molecule and check it's not empty
        smiles = self.get_smiles()
        if smiles == "":
            raise ValueError("SMILES string for this glycan is empty, check if the IUPAC is convertable.")

        # generate molecule with RDKit and check it's not None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Generated SMILES is invalid, rdkit couldn't read it in.")

        # compute the statistics and return the results inplace
        return {
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "weight": ExactMolWt(mol),
            "atoms": len(mol.GetAtoms()),
            "bonds": len(mol.GetBonds()),
            "rings": len(mol.GetRingInfo().AtomRings()),
            "monomers": len(self.parse_tree.nodes),
            "types": dict(Counter([self.parse_tree.nodes[n]["type"].get_name(True) for n in self.parse_tree.nodes])),
            "root": self.parse_tree.nodes[0]["type"].get_name(True),
            "leaves": [self.parse_tree.nodes[n]["type"].get_name(True) for n, d in self.parse_tree.out_degree() if
                       d == 0],
            "depth": max([v for k, v in nx.shortest_path_length(self.parse_tree, 0).items()]),
        }

    def count(self, glycan, match_all_fg=False, match_some_fg=False, match_edges=False, match_nodes=False,
              match_leaves=False, match_root=False):
        """Match a glycan against a query molecule and return the number of hits. This matching can be restricted by
        setting some flags introducing additional conditions of the matches.

        This matching does not include the configuration (alpha/beta/undefined) of the root monomer of the query. So
        query "Gal" will result a hit in "GalNAc6S b" but neither do "Gal a" or "Gal b".

        Args:
            glycan (Union[str, 'Glycan']): query glycan to be matched against the monomers of this glycan
            match_all_fg (bool): flag indicating to match all fgs of the query glycan to all fgs of a monomer
            match_some_fg (bool): flag indicating to match all fgs of the query glycan to some fgs of a monomer
            match_edges (bool): flag indicating to also match edges
            match_nodes (bool): flag indicating to match against all nodes
            match_leaves (bool): flag indicating to match against the leaf monomers only
            match_root (bool): flag indicating to match against the root monomer only

        Returns:
            The number of matches of the query in this glycan under the given conditions
        """
        if sum([match_nodes, match_leaves, match_root]) != 1:
            raise ValueError("Exactly one of match_nodes, match_leaves, match_root has to be True.")

        if isinstance(glycan, str):
            glycan = Glycan(glycan)

        if len(glycan.parse_tree.nodes) != 1 and (match_leaves or match_root):
            raise ValueError("Cannot match polymeric glycan against leaves of glycan. Leaves are monomers.")

        # node-match's first argument is always from first node, second argument always from second graph
        kwargs = {
            "node_match": lambda x, y: recipe_equality(x["type"], y["type"], no=True),
        }
        if match_some_fg:
            kwargs["node_match"] = lambda x, y: recipe_equality(x["type"], y["type"], some=True)
        elif match_all_fg:
            kwargs["node_match"] = lambda x, y: recipe_equality(x["type"], y["type"], every=True)
        if match_edges:
            kwargs["edge_match"] = lambda e, f: e["type"] == f["type"]

        if match_nodes:
            matcher = DiGraphMatcher(self.parse_tree, glycan.parse_tree, **kwargs)
            return len(list(matcher.subgraph_isomorphisms_iter()))
        if match_leaves:
            q = glycan.parse_tree.nodes[0]
            return sum(
                [kwargs["node_match"](self.parse_tree.nodes[n], q) for n, d in self.parse_tree.out_degree() if d == 0])
        if match_root:
            return sum([kwargs["node_match"](self.parse_tree.nodes[0], glycan.parse_tree.nodes[0])])

    def count_protonation(self, grouping):
        """
        Count the possible deprotonation sites in the final molecule.

        Args:
            grouping (bool): If True, count functional groups based on their common atom, so an SO2 group will count as
                1. Otherwise, count groups based on the protonizable oxygen atoms, so an SO2 group will count as 2.

        Returns:
            The number of possible deprotonations in the molecule.
        """
        # generate smiles for this molecule and check it's not empty
        smiles = self.get_smiles()
        if smiles == "":
            raise ValueError("SMILES string for this glycan is empty, check if the IUPAC is convertable.")

        # generate molecule with RDKit and check it's not None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Generated SMILES is invalid, rdkit couldn't read it in.")

        # iterate over deprotonatable functional groups in every form
        count = 0
        for core, group in [
            (6, [("C(O)=O", 1)]),
            (15, [("P(=O)(O)O", 3), ("P(=O)O", 2), ("P=O", 1)]),
            (16, [("S(=O)(=O)O", 2), ("S(=O)(=O)", 1)])
        ]:
            matched_atoms = set()
            for g, val in group:
                # compute the matches against this glycan
                matches = mol.GetSubstructMatches(Chem.MolFromSmiles(g))
                for match in matches:
                    for aid in match:
                        # find the core atom of each match, add it to the list of  covered atoms and increase the count
                        if mol.GetAtomWithIdx(aid).GetAtomicNum() == core and aid not in matched_atoms:
                            matched_atoms.add(aid)
                            # do some trick to either count groups to be deprotonated or possibly chargeable atoms
                            count += val ** (1 - grouping)
        return count

    def count_functional_groups(self, groups):
        """
        Count the number of the provided functional group in the final molecule.

        Args:
            groups (Union[str, List[str]]): each string has to be a valid functional group representable by this tool,
                a valid SMILES string or a valid SMARTS string

        Returns:
            The number of matches of all functional groups. This count might overlap in the matched atoms.
        """
        # in case of string input, put it into a one-element list
        if not isinstance(groups, list):
            groups = [groups]

        # generate smiles for this molecule and check it's not empty
        smiles = self.get_smiles()
        if smiles == "":
            raise ValueError("SMILES string for this glycan is empty, check if the IUPAC is convertable.")

        # generate molecule with RDKit and check it's not None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Generated SMILES is invalid, rdkit couldn't read it in.")

        fgs = []
        for group in groups:
            # convert the functional group into a RDKit molecule
            if group in functional_groups:
                tmp = Chem.MolFromSmiles(functional_groups[group])
            else:
                tmp = Chem.MolFromSmiles(group, sanitize=False)

            # check the functional groups for validity
            if tmp is None:
                raise ValueError(f"The functional group {group} cannot be parsed into a molecule. "
                                 f"Please make sure, it's a valid SMILES string or an implemented functional group.")
            fgs.append(tmp)

        # sum over matches of functional groups
        count = 0
        hydro_mol = Chem.rdmolops.AddHs(mol)
        for g in fgs:
            matched = set()
            if 1 in set(a.GetAtomicNum() for a in g.GetAtoms()):
                matches = hydro_mol.GetSubstructMatches(g)
            else:
                matches = mol.GetSubstructMatches(g)
            for match in matches:
                match = set(match)
                if len(matched.intersection(match)) == 0:
                    count += 1
                    matched = matched.union(match)
        return count

    def get_smiles(self):
        """
        Request the SMILES string of the parsed molecule.

        Returns:
            Generated SMILES string
        """
        # return an empty SMILES if the output is required to represent all modifications, but it actually wouldn't
        if not self.tree_only and self.tree_full != self.full:
            return ""

        if self.glycan_smiles is None:
            self.parse_tree, self.tree_full = TreeWalker(self.factory, False).parse(self.grammar_tree)
            self.glycan_smiles = Merger(self.factory).merge(self.parse_tree, self.root_orientation, start=self.start)
        # self.glycan_smiles = self.glycan_smiles.replace("At", "O-")
        return self.glycan_smiles

    def get_tree(self):
        """
        Request the tree parsed from the IUPAC in this instance.

        Returns:
            The parsed tree with the single monomers in the nodes.
        """
        return self.parse_tree

    def save_dot(self, output, horizontal=False):
        """
        Save the tree structure of the encoded glycan molecule into a dot file visualizing the graph of monomers.

        Args:
            output (str): path to store the DOT file in
            horizontal (bool): Show graph in horizontal orientation from left to right

        Returns:
            pydot graph object containing the graph
        """
        if horizontal:
            graph = pydot.Dot("iupac_tree", rankdir="LR")
        else:
            graph = pydot.Dot("iupac_tree")
        for node in range(len(self.parse_tree.nodes)):
            graph.add_node(pydot.Node(node, label=self.parse_tree.nodes[node]["type"].get_name(full=True)))
        for edge in self.parse_tree.edges():
            graph.add_edge(pydot.Edge(*edge[::-1], label=self.parse_tree.get_edge_data(*edge)["type"]))
        graph.write(output)
        return graph

    def create_snfg_img(self, filepath, **kwargs):
        """
        Create an image representation for a glycan using the SNFG symbols. The final image will not have a fixed size,
        but it's size adapts to the shape of the glycan. The width will be (max_depth + 1) * kwargs['width'] and the
        height will depend on the branching structure of the glycan.

        Args:
            filepath (str): path where to store the image
            **kwargs:
                - width (int): scaling factor for the width in the image generation
                - height (int): scaling factor for the width in the image generation
                - stroke (int): stroke size to be used when drawing the lines of the SNFG symbols
                - line (int): width of a line that connected two monomer-representing symbols.

        Returns:
            PIL image representation using the SNFG symbols
        """
        # set up the parameters for visualizing
        params = {
            "width": 100,
            "height": 100,
            "stroke": 2,
            "line": 3,
        }
        params.update(**kwargs)
        params["box"] = min(params["width"], params["height"])

        # generate the visualization tree representation of the glycan
        tree = Tree(self).assign_coords()
        width, height = tree.get_bounds()[2:]

        # create and draw the image
        img = Image.new(
            mode="RGB",
            size=(int((width + 1) * params["width"]), int((height + 1) * params["height"])),
            color=(255, 255, 255)
        )
        draw_img = ImageDraw.Draw(img)
        tree.draw_edges(draw_img, **params)
        tree.draw_nodes(draw_img, self, **params)

        # save the image and return is for further operations
        img.save(filepath)
        img.show()
        return img

    def __parse(self):
        """
        Adapter on the Lexer and Parser generated by ANTLR based on glyles/grammar/Glycan.g4.

        Returns:
            Nothing
        """
        # catch the prints of antlr to stderr to check if during parsing an error occurred and the glycan is invalid
        log = []

        class Writer(object):
            @staticmethod
            def write(data):
                log.append(data)

        old_err = sys.stderr
        sys.stderr = Writer()

        # parse the remaining structure description following the grammar, also add the dummy characters
        if not isinstance(self.iupac, str):
            raise ParseError("Only string input can be parsed: " + str(self.iupac))
        stream = InputStream(data='{' + self.iupac + '}')
        lexer = GlycanLexer(stream)
        token = CommonTokenStream(lexer)
        parser = GlycanParser(token)
        self.grammar_tree = parser.start()

        sys.stderr = old_err

        # if the glycan is invalid, set its structure to None and the SMILES string to empty and return
        if len(log) != 0:
            self.parse_tree = None
            self.glycan_smiles = ""
            raise ParseError("Glycan cannot be parsed:\n" + log[0])

        # walk through the AST and parse the AST into a networkx representation of the glycan.
        self.parse_tree, self.tree_full = TreeWalker(self.factory, self.tree_only).parse(self.grammar_tree)

        # if the glycan should be parsed immediately, do so
        if not self.tree_only and self.tree_full == self.full:
            self.glycan_smiles = Merger(self.factory).merge(self.parse_tree, self.root_orientation, start=self.start)


def main():
    iupac = "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-3)[Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-6)]" \
            "Gal(b1-3)[GlcNAc(a1-4)Gal(b1-4)GlcNAc6S(b1-6)]GalNAc"
    iupac2 = "IdoA2S(b1-4)Fuc2S(a1-4)[Par2S(a1-4)Xyl2S(a1-2)][All2S(a1-3)TalNAc3S(b1-6)GulN3S(a1-3)6dAltNAc3S(b1-3)]Sia6S(a1-4)Pse2S(a1-4)[Gal2S(a1-3)]Bac2S(a1-4)Tag2S"
    glycan = Glycan(iupac2, tree_only=True)
    glycan.create_snfg_img("test.png", **{
        "width": 300,
        "height": 100,
        "stroke": 2,
        "line": 3,
    })


if __name__ == '__main__':
    main()
