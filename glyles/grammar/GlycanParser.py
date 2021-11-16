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
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\t")
        buf.write("\'\4\2\t\2\4\3\t\3\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2")
        buf.write("\3\2\3\2\3\2\5\2\23\n\2\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3")
        buf.write("\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\5\3%\n\3\3\3\2\2\4")
        buf.write("\2\4\2\2\2*\2\22\3\2\2\2\4$\3\2\2\2\6\7\5\4\3\2\7\b\7")
        buf.write("\6\2\2\b\t\7\3\2\2\t\n\7\b\2\2\n\23\3\2\2\2\13\f\5\4\3")
        buf.write("\2\f\r\7\6\2\2\r\23\3\2\2\2\16\17\7\6\2\2\17\20\7\3\2")
        buf.write("\2\20\23\7\b\2\2\21\23\7\6\2\2\22\6\3\2\2\2\22\13\3\2")
        buf.write("\2\2\22\16\3\2\2\2\22\21\3\2\2\2\23\3\3\2\2\2\24\25\7")
        buf.write("\6\2\2\25%\7\7\2\2\26\27\7\6\2\2\27\30\7\7\2\2\30%\5\4")
        buf.write("\3\2\31\32\7\4\2\2\32\33\5\4\3\2\33\34\7\5\2\2\34%\3\2")
        buf.write("\2\2\35\36\7\6\2\2\36\37\7\7\2\2\37 \7\4\2\2 !\5\4\3\2")
        buf.write("!\"\7\5\2\2\"#\5\4\3\2#%\3\2\2\2$\24\3\2\2\2$\26\3\2\2")
        buf.write("\2$\31\3\2\2\2$\35\3\2\2\2%\5\3\2\2\2\4\22$")
        return buf.getvalue()


class GlycanParser ( Parser ):

    grammarFileName = "Glycan.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "' '", "'['", "']'" ]

    symbolicNames = [ "<INVALID>", "<INVALID>", "<INVALID>", "<INVALID>", 
                      "SAC", "CON", "TYPE", "NUM" ]

    RULE_start = 0
    RULE_branch = 1

    ruleNames =  [ "start", "branch" ]

    EOF = Token.EOF
    T__0=1
    T__1=2
    T__2=3
    SAC=4
    CON=5
    TYPE=6
    NUM=7

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

        def TYPE(self):
            return self.getToken(GlycanParser.TYPE, 0)

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
            self.state = 16
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,0,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 4
                self.branch()
                self.state = 5
                self.match(GlycanParser.SAC)
                self.state = 6
                self.match(GlycanParser.T__0)
                self.state = 7
                self.match(GlycanParser.TYPE)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 9
                self.branch()
                self.state = 10
                self.match(GlycanParser.SAC)
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 12
                self.match(GlycanParser.SAC)
                self.state = 13
                self.match(GlycanParser.T__0)
                self.state = 14
                self.match(GlycanParser.TYPE)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 15
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

        def SAC(self):
            return self.getToken(GlycanParser.SAC, 0)

        def CON(self):
            return self.getToken(GlycanParser.CON, 0)

        def branch(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GlycanParser.BranchContext)
            else:
                return self.getTypedRuleContext(GlycanParser.BranchContext,i)


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
            self.state = 34
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,1,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 18
                self.match(GlycanParser.SAC)
                self.state = 19
                self.match(GlycanParser.CON)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 20
                self.match(GlycanParser.SAC)
                self.state = 21
                self.match(GlycanParser.CON)
                self.state = 22
                self.branch()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 23
                self.match(GlycanParser.T__1)
                self.state = 24
                self.branch()
                self.state = 25
                self.match(GlycanParser.T__2)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 27
                self.match(GlycanParser.SAC)
                self.state = 28
                self.match(GlycanParser.CON)
                self.state = 29
                self.match(GlycanParser.T__1)
                self.state = 30
                self.branch()
                self.state = 31
                self.match(GlycanParser.T__2)
                self.state = 32
                self.branch()
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





