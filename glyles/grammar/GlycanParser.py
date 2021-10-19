# Generated from Glycan.g4 by ANTLR 4.9
# encoding: utf-8
from antlr4 import *
from io import StringIO
import sys
if sys.version_info[1] > 5:
	from typing import TextIO
else:
	from typing.io import TextIO


def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\6")
        buf.write("\34\4\2\t\2\4\3\t\3\3\2\3\2\3\3\3\3\3\3\3\3\3\3\3\3\3")
        buf.write("\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\5\3\32\n\3")
        buf.write("\3\3\2\2\4\2\4\2\2\2\34\2\6\3\2\2\2\4\31\3\2\2\2\6\7\5")
        buf.write("\4\3\2\7\3\3\2\2\2\b\t\7\5\2\2\t\32\7\6\2\2\n\13\7\5\2")
        buf.write("\2\13\f\7\6\2\2\f\32\5\4\3\2\r\16\7\3\2\2\16\17\5\4\3")
        buf.write("\2\17\20\7\4\2\2\20\32\3\2\2\2\21\22\7\5\2\2\22\23\7\6")
        buf.write("\2\2\23\24\7\3\2\2\24\25\5\4\3\2\25\26\7\4\2\2\26\27\7")
        buf.write("\5\2\2\27\30\7\6\2\2\30\32\3\2\2\2\31\b\3\2\2\2\31\n\3")
        buf.write("\2\2\2\31\r\3\2\2\2\31\21\3\2\2\2\32\5\3\2\2\2\3\31")
        return buf.getvalue()


class GlycanParser ( Parser ):

    grammarFileName = "Glycan.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'['", "']'" ]

    symbolicNames = [ "<INVALID>", "<INVALID>", "<INVALID>", "SAC", "CON" ]

    RULE_start = 0
    RULE_branch = 1

    ruleNames =  [ "start", "branch" ]

    EOF = Token.EOF
    T__0=1
    T__1=2
    SAC=3
    CON=4

    def __init__(self, input:TokenStream, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.9")
        self._interp = ParserATNSimulator(self, self.atn, self.decisionsToDFA, self.sharedContextCache)
        self._predicates = None




    class StartContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def branch(self):
            return self.getTypedRuleContext(GlycanParser.BranchContext,0)


        def getRuleIndex(self):
            return GlycanParser.RULE_start

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterStart" ):
                listener.enterStart(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitStart" ):
                listener.exitStart(self)




    def start(self):

        localctx = GlycanParser.StartContext(self, self._ctx, self.state)
        self.enterRule(localctx, 0, self.RULE_start)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 4
            self.branch()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class BranchContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def SAC(self, i:int=None):
            if i is None:
                return self.getTokens(GlycanParser.SAC)
            else:
                return self.getToken(GlycanParser.SAC, i)

        def CON(self, i:int=None):
            if i is None:
                return self.getTokens(GlycanParser.CON)
            else:
                return self.getToken(GlycanParser.CON, i)

        def branch(self):
            return self.getTypedRuleContext(GlycanParser.BranchContext,0)


        def getRuleIndex(self):
            return GlycanParser.RULE_branch

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterBranch" ):
                listener.enterBranch(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitBranch" ):
                listener.exitBranch(self)




    def branch(self):

        localctx = GlycanParser.BranchContext(self, self._ctx, self.state)
        self.enterRule(localctx, 2, self.RULE_branch)
        try:
            self.state = 23
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,0,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 6
                self.match(GlycanParser.SAC)
                self.state = 7
                self.match(GlycanParser.CON)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 8
                self.match(GlycanParser.SAC)
                self.state = 9
                self.match(GlycanParser.CON)
                self.state = 10
                self.branch()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 11
                self.match(GlycanParser.T__0)
                self.state = 12
                self.branch()
                self.state = 13
                self.match(GlycanParser.T__1)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 15
                self.match(GlycanParser.SAC)
                self.state = 16
                self.match(GlycanParser.CON)
                self.state = 17
                self.match(GlycanParser.T__0)
                self.state = 18
                self.branch()
                self.state = 19
                self.match(GlycanParser.T__1)
                self.state = 20
                self.match(GlycanParser.SAC)
                self.state = 21
                self.match(GlycanParser.CON)
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





