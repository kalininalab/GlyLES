# Generated from Glycan.g4 by ANTLR 4.9
from antlr4 import *
if __name__ is not None and "." in __name__:
    from .GlycanParser import GlycanParser
else:
    from GlycanParser import GlycanParser

# This class defines a complete listener for a parse tree produced by GlycanParser.
class GlycanListener(ParseTreeListener):

    # Enter a parse tree produced by GlycanParser#start.
    def enterStart(self, ctx:GlycanParser.StartContext):
        pass

    # Exit a parse tree produced by GlycanParser#start.
    def exitStart(self, ctx:GlycanParser.StartContext):
        pass


    # Enter a parse tree produced by GlycanParser#branch.
    def enterBranch(self, ctx:GlycanParser.BranchContext):
        pass

    # Exit a parse tree produced by GlycanParser#branch.
    def exitBranch(self, ctx:GlycanParser.BranchContext):
        pass


    # Enter a parse tree produced by GlycanParser#glycan.
    def enterGlycan(self, ctx:GlycanParser.GlycanContext):
        pass

    # Exit a parse tree produced by GlycanParser#glycan.
    def exitGlycan(self, ctx:GlycanParser.GlycanContext):
        pass


    # Enter a parse tree produced by GlycanParser#deriv.
    def enterDeriv(self, ctx:GlycanParser.DerivContext):
        pass

    # Exit a parse tree produced by GlycanParser#deriv.
    def exitDeriv(self, ctx:GlycanParser.DerivContext):
        pass



del GlycanParser