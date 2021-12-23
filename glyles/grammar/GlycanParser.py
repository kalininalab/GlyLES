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
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\13")
        buf.write("/\4\2\t\2\4\3\t\3\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2")
        buf.write("\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\5\2\33\n")
        buf.write("\2\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3")
        buf.write("\3\3\3\3\3\3\5\3-\n\3\3\3\2\2\4\2\4\2\2\2\62\2\32\3\2")
        buf.write("\2\2\4,\3\2\2\2\6\7\7\3\2\2\7\b\5\4\3\2\b\t\7\b\2\2\t")
        buf.write("\n\7\4\2\2\n\13\7\n\2\2\13\f\7\5\2\2\f\33\3\2\2\2\r\16")
        buf.write("\7\3\2\2\16\17\5\4\3\2\17\20\7\b\2\2\20\21\7\5\2\2\21")
        buf.write("\33\3\2\2\2\22\23\7\3\2\2\23\24\7\b\2\2\24\25\7\4\2\2")
        buf.write("\25\26\7\n\2\2\26\33\7\5\2\2\27\30\7\3\2\2\30\31\7\b\2")
        buf.write("\2\31\33\7\5\2\2\32\6\3\2\2\2\32\r\3\2\2\2\32\22\3\2\2")
        buf.write("\2\32\27\3\2\2\2\33\3\3\2\2\2\34\35\7\b\2\2\35-\7\t\2")
        buf.write("\2\36\37\7\b\2\2\37 \7\t\2\2 -\5\4\3\2!\"\7\6\2\2\"#\5")
        buf.write("\4\3\2#$\7\7\2\2$-\3\2\2\2%&\7\b\2\2&\'\7\t\2\2\'(\7\6")
        buf.write("\2\2()\5\4\3\2)*\7\7\2\2*+\5\4\3\2+-\3\2\2\2,\34\3\2\2")
        buf.write("\2,\36\3\2\2\2,!\3\2\2\2,%\3\2\2\2-\5\3\2\2\2\4\32,")
        return buf.getvalue()


class GlycanParser ( Parser ):

    grammarFileName = "Glycan.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'{'", "' '", "'}'", "'['", "']'" ]

    symbolicNames = [ "<INVALID>", "<INVALID>", "<INVALID>", "<INVALID>", 
                      "<INVALID>", "<INVALID>", "SAC", "CON", "TYPE", "NUM" ]

    RULE_start = 0
    RULE_branch = 1

    ruleNames =  [ "start", "branch" ]

    EOF = Token.EOF
    T__0=1
    T__1=2
    T__2=3
    T__3=4
    T__4=5
    SAC=6
    CON=7
    TYPE=8
    NUM=9

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
            self.state = 24
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,0,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 4
                self.match(GlycanParser.T__0)
                self.state = 5
                self.branch()
                self.state = 6
                self.match(GlycanParser.SAC)
                self.state = 7
                self.match(GlycanParser.T__1)
                self.state = 8
                self.match(GlycanParser.TYPE)
                self.state = 9
                self.match(GlycanParser.T__2)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 11
                self.match(GlycanParser.T__0)
                self.state = 12
                self.branch()
                self.state = 13
                self.match(GlycanParser.SAC)
                self.state = 14
                self.match(GlycanParser.T__2)
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 16
                self.match(GlycanParser.T__0)
                self.state = 17
                self.match(GlycanParser.SAC)
                self.state = 18
                self.match(GlycanParser.T__1)
                self.state = 19
                self.match(GlycanParser.TYPE)
                self.state = 20
                self.match(GlycanParser.T__2)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 21
                self.match(GlycanParser.T__0)
                self.state = 22
                self.match(GlycanParser.SAC)
                self.state = 23
                self.match(GlycanParser.T__2)
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
            self.state = 42
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,1,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 26
                self.match(GlycanParser.SAC)
                self.state = 27
                self.match(GlycanParser.CON)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 28
                self.match(GlycanParser.SAC)
                self.state = 29
                self.match(GlycanParser.CON)
                self.state = 30
                self.branch()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 31
                self.match(GlycanParser.T__3)
                self.state = 32
                self.branch()
                self.state = 33
                self.match(GlycanParser.T__4)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 35
                self.match(GlycanParser.SAC)
                self.state = 36
                self.match(GlycanParser.CON)
                self.state = 37
                self.match(GlycanParser.T__3)
                self.state = 38
                self.branch()
                self.state = 39
                self.match(GlycanParser.T__4)
                self.state = 40
                self.branch()
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





