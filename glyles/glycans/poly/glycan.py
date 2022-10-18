import sys
from typing import Union, List

from networkx.algorithms.isomorphism import DiGraphMatcher
import pydot
from antlr4 import InputStream, CommonTokenStream
from rdkit import Chem

from glyles.glycans.factory.factory import MonomerFactory
from glyles.glycans.mono.reactor import functional_groups
from glyles.glycans.poly.merger import Merger
from glyles.glycans.poly.walker import TreeWalker
from glyles.glycans.utils import ParseError, find_isomorphism_nx
from glyles.grammar.GlycanLexer import GlycanLexer
from glyles.grammar.GlycanParser import GlycanParser


def compare_smiles(c, s):
    Chem.Kekulize(c)
    Chem.Kekulize(s)

    ssmiles = Chem.MolToSmiles(s, kekuleSmiles=True)
    csmiles = Chem.MolToSmiles(c, kekuleSmiles=True)
    return csmiles == ssmiles


def recipe_equality(glycan, query, no=False, some=False, every=False):
    """
    Checking if two monomers are considered isomorphic under the given mode of isomorphism.
    Note:
        node-match's first argument is always from first node, second argument always from second graph

    Args:
        glycan (Monomer):
        query (Monomer):
        no (bool):
        some (bool):
        every (bool):

    Returns:

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
        raise NotImplementedError("Matching to a subset of functional groups is not yet implemented. Coming soon.")
        # if a query's fgs have to be a subset of glycan's fgs, do some magic
        recipe_glycan, recipe_query = glycan.get_recipe(), query.get_recipe()
        if recipe_glycan[list(zip(*recipe_glycan))[1].index(GlycanLexer.SAC)] != \
                recipe_query[list(zip(*recipe_query))[1].index(GlycanLexer.SAC)]:
            return False
        iso = find_isomorphism_nx(glycan, query, query.name, query.c1_find)
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
            tree_only (bool): Flag indicating to only parse the tree of glycans and not the modifications
            full (bool): Flag indicating that only fully convertible glycans should be returned, i.e. all modifications
                such as 3-Anhydro-[...] are also present in the SMILES
        """
        self.iupac = iupac
        self.parse_tree = None
        self.glycan_smiles = None
        self.root_orientation = root_orientation
        self.start = start
        self.tree_only = tree_only
        self.factory = MonomerFactory()
        self.full = full
        self.tree_full = True
        self.__parse()

    def count(
            self,
            glycan,
            match_all_fg=False,
            match_some_fg=False,
            match_edges=False,
            match_nodes=False,
            match_leaves=False,
            match_root=False,
    ):
        """
        Match a glycan against a query molecule and return the number of hits. This matching can be restricted by 
        setting some flags introducing additional conditions of the matches.

        This matching does not include the configuration (alpha/beta/undefined) of the root monomer of the query. So
        query "Gal" will result a hit in "GalNAc6S b" but neither do "Gal a" or "Gal b".

        Args:
            glycan (Union[str, Glycan]): query glycan to be matched against the monomers of this glycan
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
            glycan = Glycan(glycan, full=False)

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
            return sum([kwargs["node_match"](self.parse_tree.nodes[n], q) for n, d in self.parse_tree.out_degree() if d == 0])
        if match_root:
            return sum([kwargs["node_match"](self.parse_tree.nodes[0], glycan.parse_tree.nodes[0])])

    def count_protonation(self, groups):
        """
        Count the possible deprotonation sites in the final molecule.

        Args:
            groups (bool): If True, count functional groups that can be deprotonated; otherwise, count possible
                deprotonations

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
                            # do some trick to either count groups to be deprotonated or possibly chargable atoms
                            count += val ** groups
        return count

    def count_functional_groups(self, groups: Union[str, List[str]]):
        """
        Count the number of the provided functional group in the final molecule.

        Args:
            groups: string of a specific group to find or a list of strings.
                Those strings have to be valid SMILES strings.

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
                tmp = Chem.MolFromSmiles(group)
            # check the functional groups for validity
            if tmp is None:
                raise ValueError(f"The functional group {group} cannot be parsed into a molecule. "
                                 f"Please make sure, it's a valid SMILES string or an implemented functional group.")
            fgs.append(tmp)

        # sum over matches of functional groups
        return sum([len(mol.GetSubstructMatches(g)) for g in fgs])

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
            self.glycan_smiles = Merger(self.factory).merge(self.parse_tree, self.root_orientation, start=self.start)
        self.glycan_smiles = self.glycan_smiles.replace("At", "O-")
        return self.glycan_smiles

    def get_tree(self):
        """
        Request the tree parsed from the IUPAC in this instance.

        Returns:
            The parsed tree with the single monomers in the nodes
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
        tree = parser.start()

        sys.stderr = old_err

        # if the glycan is invalid, set its structure to None and the SMILES string to empty and return
        if len(log) != 0:
            self.parse_tree = None
            self.glycan_smiles = ""
            raise ParseError("Glycan cannot be parsed:\n" + log[0])

        # walk through the AST and parse the AST into a networkx representation of the glycan.
        self.parse_tree, self.tree_full = TreeWalker(self.factory, self.tree_only).parse(tree)

        # if the glycan should be parsed immediately, do so
        if not self.tree_only and self.tree_full == self.full:
            self.glycan_smiles = Merger(self.factory).merge(self.parse_tree, self.root_orientation, start=self.start)