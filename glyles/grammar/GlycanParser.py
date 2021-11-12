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
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\b")
        buf.write(" \4\2\t\2\4\3\t\3\3\2\3\2\3\2\3\2\5\2\13\n\2\3\3\3\3\3")
        buf.write("\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3")
        buf.write("\3\3\5\3\36\n\3\3\3\2\2\4\2\4\2\2\2!\2\n\3\2\2\2\4\35")
        buf.write("\3\2\2\2\6\7\5\4\3\2\7\b\7\5\2\2\b\13\3\2\2\2\t\13\7\5")
        buf.write("\2\2\n\6\3\2\2\2\n\t\3\2\2\2\13\3\3\2\2\2\f\r\7\5\2\2")
        buf.write("\r\36\7\6\2\2\16\17\7\5\2\2\17\20\7\6\2\2\20\36\5\4\3")
        buf.write("\2\21\22\7\3\2\2\22\23\5\4\3\2\23\24\7\4\2\2\24\36\3\2")
        buf.write("\2\2\25\26\7\5\2\2\26\27\7\6\2\2\27\30\7\3\2\2\30\31\5")
        buf.write("\4\3\2\31\32\7\4\2\2\32\33\7\5\2\2\33\34\7\6\2\2\34\36")
        buf.write("\3\2\2\2\35\f\3\2\2\2\35\16\3\2\2\2\35\21\3\2\2\2\35\25")
        buf.write("\3\2\2\2\36\5\3\2\2\2\4\n\35")
        return buf.getvalue()


class GlycanParser ( Parser ):

    grammarFileName = "Glycan.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'['", "']'" ]

    symbolicNames = [ "<INVALID>", "<INVALID>", "<INVALID>", "SAC", "CON", 
                      "TYPE", "NUM" ]

    RULE_start = 0
    RULE_branch = 1

    ruleNames =  [ "start", "branch" ]

    EOF = Token.EOF
    T__0=1
    T__1=2
    SAC=3
    CON=4
    TYPE=5
    NUM=6

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


        def SAC(self):
            return self.getToken(GlycanParser.SAC, 0)

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
            self.state = 8
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,0,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 4
                self.branch()
                self.state = 5
                self.match(GlycanParser.SAC)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 7
                self.match(GlycanParser.SAC)
                pass


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
            self.state = 27
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,1,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 10
                self.match(GlycanParser.SAC)
                self.state = 11
                self.match(GlycanParser.CON)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 12
                self.match(GlycanParser.SAC)
                self.state = 13
                self.match(GlycanParser.CON)
                self.state = 14
                self.branch()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 15
                self.match(GlycanParser.T__0)
                self.state = 16
                self.branch()
                self.state = 17
                self.match(GlycanParser.T__1)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 19
                self.match(GlycanParser.SAC)
                self.state = 20
                self.match(GlycanParser.CON)
                self.state = 21
                self.match(GlycanParser.T__0)
                self.state = 22
                self.branch()
                self.state = 23
                self.match(GlycanParser.T__1)
                self.state = 24
                self.match(GlycanParser.SAC)
                self.state = 25
                self.match(GlycanParser.CON)
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





