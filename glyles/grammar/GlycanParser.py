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
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\f")
        buf.write("?\4\2\t\2\4\3\t\3\4\4\t\4\3\2\3\2\3\2\3\2\3\2\3\2\3\2")
        buf.write("\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3\2\3")
        buf.write("\2\3\2\5\2\37\n\2\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3")
        buf.write("\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\5\3\63\n\3\3\4\3")
        buf.write("\4\3\4\3\4\3\4\3\4\3\4\3\4\5\4=\n\4\3\4\2\2\5\2\4\6\2")
        buf.write("\2\2D\2\36\3\2\2\2\4\62\3\2\2\2\6<\3\2\2\2\b\t\7\3\2\2")
        buf.write("\t\n\5\4\3\2\n\13\5\6\4\2\13\f\7\4\2\2\f\r\7\n\2\2\r\16")
        buf.write("\7\5\2\2\16\37\3\2\2\2\17\20\7\3\2\2\20\21\5\4\3\2\21")
        buf.write("\22\5\6\4\2\22\23\7\5\2\2\23\37\3\2\2\2\24\25\7\3\2\2")
        buf.write("\25\26\5\6\4\2\26\27\7\4\2\2\27\30\7\n\2\2\30\31\7\5\2")
        buf.write("\2\31\37\3\2\2\2\32\33\7\3\2\2\33\34\5\6\4\2\34\35\7\5")
        buf.write("\2\2\35\37\3\2\2\2\36\b\3\2\2\2\36\17\3\2\2\2\36\24\3")
        buf.write("\2\2\2\36\32\3\2\2\2\37\3\3\2\2\2 !\5\6\4\2!\"\7\t\2\2")
        buf.write("\"\63\3\2\2\2#$\5\6\4\2$%\7\t\2\2%&\5\4\3\2&\63\3\2\2")
        buf.write("\2\'(\7\6\2\2()\5\4\3\2)*\7\7\2\2*\63\3\2\2\2+,\5\6\4")
        buf.write("\2,-\7\t\2\2-.\7\6\2\2./\5\4\3\2/\60\7\7\2\2\60\61\5\4")
        buf.write("\3\2\61\63\3\2\2\2\62 \3\2\2\2\62#\3\2\2\2\62\'\3\2\2")
        buf.write("\2\62+\3\2\2\2\63\5\3\2\2\2\64\65\7\b\2\2\65\66\7\13\2")
        buf.write("\2\66=\7\n\2\2\678\7\b\2\28=\7\n\2\29:\7\b\2\2:=\7\13")
        buf.write("\2\2;=\7\b\2\2<\64\3\2\2\2<\67\3\2\2\2<9\3\2\2\2<;\3\2")
        buf.write("\2\2=\7\3\2\2\2\5\36\62<")
        return buf.getvalue()


class GlycanParser ( Parser ):

    grammarFileName = "Glycan.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'{'", "' '", "'}'", "'['", "']'" ]

    symbolicNames = [ "<INVALID>", "<INVALID>", "<INVALID>", "<INVALID>", 
                      "<INVALID>", "<INVALID>", "SAC", "CON", "TYPE", "RING", 
                      "NUM" ]

    RULE_start = 0
    RULE_branch = 1
    RULE_glycan = 2

    ruleNames =  [ "start", "branch", "glycan" ]

    EOF = Token.EOF
    T__0=1
    T__1=2
    T__2=3
    T__3=4
    T__4=5
    SAC=6
    CON=7
    TYPE=8
    RING=9
    NUM=10

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


        def glycan(self):
            return self.getTypedRuleContext(GlycanParser.GlycanContext,0)


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
            self.state = 28
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,0,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 6
                self.match(GlycanParser.T__0)
                self.state = 7
                self.branch()
                self.state = 8
                self.glycan()
                self.state = 9
                self.match(GlycanParser.T__1)
                self.state = 10
                self.match(GlycanParser.TYPE)
                self.state = 11
                self.match(GlycanParser.T__2)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 13
                self.match(GlycanParser.T__0)
                self.state = 14
                self.branch()
                self.state = 15
                self.glycan()
                self.state = 16
                self.match(GlycanParser.T__2)
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 18
                self.match(GlycanParser.T__0)
                self.state = 19
                self.glycan()
                self.state = 20
                self.match(GlycanParser.T__1)
                self.state = 21
                self.match(GlycanParser.TYPE)
                self.state = 22
                self.match(GlycanParser.T__2)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 24
                self.match(GlycanParser.T__0)
                self.state = 25
                self.glycan()
                self.state = 26
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

        def glycan(self):
            return self.getTypedRuleContext(GlycanParser.GlycanContext,0)


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
            self.state = 48
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,1,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 30
                self.glycan()
                self.state = 31
                self.match(GlycanParser.CON)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 33
                self.glycan()
                self.state = 34
                self.match(GlycanParser.CON)
                self.state = 35
                self.branch()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 37
                self.match(GlycanParser.T__3)
                self.state = 38
                self.branch()
                self.state = 39
                self.match(GlycanParser.T__4)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 41
                self.glycan()
                self.state = 42
                self.match(GlycanParser.CON)
                self.state = 43
                self.match(GlycanParser.T__3)
                self.state = 44
                self.branch()
                self.state = 45
                self.match(GlycanParser.T__4)
                self.state = 46
                self.branch()
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class GlycanContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def SAC(self):
            return self.getToken(GlycanParser.SAC, 0)

        def RING(self):
            return self.getToken(GlycanParser.RING, 0)

        def TYPE(self):
            return self.getToken(GlycanParser.TYPE, 0)

        def getRuleIndex(self):
            return GlycanParser.RULE_glycan

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterGlycan" ):
                listener.enterGlycan(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitGlycan" ):
                listener.exitGlycan(self)




    def glycan(self):

        localctx = GlycanParser.GlycanContext(self, self._ctx, self.state)
        self.enterRule(localctx, 4, self.RULE_glycan)
        try:
            self.state = 58
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,2,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 50
                self.match(GlycanParser.SAC)
                self.state = 51
                self.match(GlycanParser.RING)
                self.state = 52
                self.match(GlycanParser.TYPE)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 53
                self.match(GlycanParser.SAC)
                self.state = 54
                self.match(GlycanParser.TYPE)
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 55
                self.match(GlycanParser.SAC)
                self.state = 56
                self.match(GlycanParser.RING)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 57
                self.match(GlycanParser.SAC)
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





