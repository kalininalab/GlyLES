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
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\n")
        buf.write("&\4\2\t\2\4\3\t\3\4\4\t\4\3\2\3\2\3\3\3\3\3\3\3\3\3\3")
        buf.write("\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\5")
        buf.write("\3\35\n\3\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\2\2\5\2\4\6")
        buf.write("\2\2\2%\2\b\3\2\2\2\4\34\3\2\2\2\6\36\3\2\2\2\b\t\5\4")
        buf.write("\3\2\t\3\3\2\2\2\n\13\7\b\2\2\13\35\5\6\4\2\f\r\7\b\2")
        buf.write("\2\r\16\5\6\4\2\16\17\5\4\3\2\17\35\3\2\2\2\20\21\7\3")
        buf.write("\2\2\21\22\5\4\3\2\22\23\7\4\2\2\23\35\3\2\2\2\24\25\7")
        buf.write("\b\2\2\25\26\5\6\4\2\26\27\7\3\2\2\27\30\5\4\3\2\30\31")
        buf.write("\7\4\2\2\31\32\7\b\2\2\32\33\5\6\4\2\33\35\3\2\2\2\34")
        buf.write("\n\3\2\2\2\34\f\3\2\2\2\34\20\3\2\2\2\34\24\3\2\2\2\35")
        buf.write("\5\3\2\2\2\36\37\7\5\2\2\37 \7\t\2\2 !\7\n\2\2!\"\7\6")
        buf.write("\2\2\"#\7\n\2\2#$\7\7\2\2$\7\3\2\2\2\3\34")
        return buf.getvalue()


class GlycanParser ( Parser ):

    grammarFileName = "Glycan.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'['", "']'", "'('", "'-'", "')'" ]

    symbolicNames = [ "<INVALID>", "<INVALID>", "<INVALID>", "<INVALID>", 
                      "<INVALID>", "<INVALID>", "SAC", "TYPE", "NUM" ]

    RULE_start = 0
    RULE_branch = 1
    RULE_con = 2

    ruleNames =  [ "start", "branch", "con" ]

    EOF = Token.EOF
    T__0=1
    T__1=2
    T__2=3
    T__3=4
    T__4=5
    SAC=6
    TYPE=7
    NUM=8

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
            self.state = 6
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

        def con(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GlycanParser.ConContext)
            else:
                return self.getTypedRuleContext(GlycanParser.ConContext,i)


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
            self.state = 26
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,0,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 8
                self.match(GlycanParser.SAC)
                self.state = 9
                self.con()
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 10
                self.match(GlycanParser.SAC)
                self.state = 11
                self.con()
                self.state = 12
                self.branch()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 14
                self.match(GlycanParser.T__0)
                self.state = 15
                self.branch()
                self.state = 16
                self.match(GlycanParser.T__1)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 18
                self.match(GlycanParser.SAC)
                self.state = 19
                self.con()
                self.state = 20
                self.match(GlycanParser.T__0)
                self.state = 21
                self.branch()
                self.state = 22
                self.match(GlycanParser.T__1)
                self.state = 23
                self.match(GlycanParser.SAC)
                self.state = 24
                self.con()
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class ConContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def TYPE(self):
            return self.getToken(GlycanParser.TYPE, 0)

        def NUM(self, i:int=None):
            if i is None:
                return self.getTokens(GlycanParser.NUM)
            else:
                return self.getToken(GlycanParser.NUM, i)

        def getRuleIndex(self):
            return GlycanParser.RULE_con

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterCon" ):
                listener.enterCon(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitCon" ):
                listener.exitCon(self)

        def __str__(self):
            return "".join(c.symbol.text for c in self.children)


    def con(self):

        localctx = GlycanParser.ConContext(self, self._ctx, self.state)
        self.enterRule(localctx, 4, self.RULE_con)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 28
            self.match(GlycanParser.T__2)
            self.state = 29
            self.match(GlycanParser.TYPE)
            self.state = 30
            self.match(GlycanParser.NUM)
            self.state = 31
            self.match(GlycanParser.T__3)
            self.state = 32
            self.match(GlycanParser.NUM)
            self.state = 33
            self.match(GlycanParser.T__4)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





