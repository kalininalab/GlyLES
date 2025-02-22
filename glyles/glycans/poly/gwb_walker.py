import networkx as nx
from antlr4 import ErrorNode, TerminalNode

from glyles.glycans.poly.walker import TreeWalker
from glyles.glycans.utils import UnreachableError
from glyles.iupac.IUPACLexer import IUPACLexer
from glyles.gwb.GWBLexer import GWBLexer
from glyles.gwb.GWBParser import GWBParser


class GWBTreeWalker(TreeWalker):
    def __init__(self, factory, tree_only):
        super().__init__(factory, tree_only)
    
    def parse(self, t):
        floating = False
        children = list(t.getChildren())
        for c in filter(lambda x: not isinstance(x, (TerminalNode, ErrorNode)), children):
            if isinstance(c, GWBParser.BranchContext):
                self.walk(c, -1, floating=floating)
                floating = True
        
        if isinstance(children[0], GWBParser.BosContext) and children[0].children[0].symbol.type == GWBLexer.FREEEND:
            self.update_recipe(0, [("-ol", GWBLexer.MOD)])
        
        for node in self.g.nodes:
            if "update" in self.g.nodes[node]:
                self.update_node(node)
        
        return self.g, self.full and len(list(nx.connected_components(self.g.to_undirected()))) == 1
    
    def walk(self, t, parent, floating: bool = False):
        """
        Args:
            t (GWBParser.BranchContext): The branch to parse.
        """
        if isinstance(t, ErrorNode):
            raise UnreachableError("ErrorNodes in parsing!")
        elif isinstance(t, TerminalNode):
            raise UnreachableError("TerminalNodes in parsing!")
        
        children = list(t.getChildren())
        if len(children) == 2:  # con deriv
            node_id = self.add_node(children[1], parent=parent)
            self.full &= self.add_edge(parent, node_id, children[0], floating)
            return node_id
        if len(children) == 4:  # ( branch ) branch
            self.walk(children[1], parent, floating)
            return self.walk(children[3], parent, floating)
        elif len(children) in {3, 6}:  # con deriv (LPAR branch RPAR)? branch
            node_id = self.add_node(children[1], parent=parent)
            if len(children) > 3:
                self.walk(children[3], node_id, floating)  # link branches to node_id, which is the upstream monosaccharide
            self.walk(children[-1], node_id, floating)
            self.full &= self.add_edge(parent, node_id, children[0], floating)
            return node_id
        raise UnreachableError("Unexpected number of children in branch!")

    def add_edge(self, parent, child, con, floating: bool = False) -> bool:
        cs = [self.context2str(c) for c in list(con.getChildren())[1:]]
        if len(cs) < 3:
            cs += ["?" for _ in range(len(cs), 3)]
        con = f"({cs[1]}{cs[2]}-{cs[0]})"
        if parent == -1 and floating:
            self.g.nodes[child]["float_edge"] = con
            return False
        if parent != -1:
            if parent == child:  # special case for branches with sulfate modification or similar
                if cs[0] != "?":
                    self.g.nodes[parent]["update"] = [(cs[0], IUPACLexer.NUM)] + self.g.nodes[parent]["update"]
                return True
            self.g.add_edge(parent, child, type=con)
            return "?" not in con and "/" not in con
        return True  # parent == -1 and not floating

    def update_node(self, node_id):
        recipe = self.g.nodes[node_id]["update"]
        num = None
        new_recipe = []
        for tok, type_ in recipe:
            if type_ == GWBLexer.NUM:
                num = tok
            else:
                if num is not None:
                    tok = num + tok
                    num = None
                new_recipe.append((tok, type_))
        p = self.g.nodes[node_id]["type"]
        monomer, full = self.factory.create(p.get_recipe() + new_recipe, p.get_config().to_string(), tree_only=self.tree_only)
        self.g.nodes[node_id]["type"] = monomer
        self.full &= full
        del self.g.nodes[node_id]["update"]
    
    def update_recipe(self, parent, recipe):
        if "update" in self.g.nodes[parent]:
            self.g.nodes[parent]["update"] += recipe
        else:
            self.g.nodes[parent]["update"] = recipe
    
    def add_node(self, node, parent, config=""):
        """
        Add a new node to the network based on the name of the represented glycan.

        Args:
            node (IUPACParser.GlycanContext): Node of the parsed tree to be parsed into a monomer
            config (str): configuration to be applied to the monomer

        Returns:
            ID of the newly added node
        """
        # get the id for the node and increase the id for the next node
        node_id = self.node_id
        self.node_id += 1

        # add the node to the network and store the enum of the glycan as attribute
        recipe = self.build_recipe(node)
        if len(set.intersection({IUPACLexer.SAC, IUPACLexer.COUNT}, set(x for _, x in recipe))) == 0:
            self.update_recipe(parent, recipe)
            return parent
        elif len(recipe) > 2 and recipe[-2][0] == ",":
            recipe = recipe[:-2] + [recipe[-1]]
        monomer, full = self.factory.create(recipe, config, tree_only=self.tree_only)
        self.g.add_node(
            node_id,
            type=monomer,
        )
        self.full &= full

        return node_id
