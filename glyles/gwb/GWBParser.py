# Generated from GWB.g4 by ANTLR 4.13.2
# encoding: utf-8
from antlr4 import *
from io import StringIO
import sys
if sys.version_info[1] > 5:
	from typing import TextIO
else:
	from typing.io import TextIO

def serializedATN():
    return [
        4,1,37,214,2,0,7,0,2,1,7,1,2,2,7,2,2,3,7,3,2,4,7,4,2,5,7,5,2,6,7,
        6,2,7,7,7,2,8,7,8,2,9,7,9,2,10,7,10,2,11,7,11,2,12,7,12,2,13,7,13,
        1,0,1,0,1,0,1,0,3,0,33,8,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,51,8,1,1,1,1,1,3,1,55,8,1,1,2,1,2,
        3,2,59,8,2,1,2,5,2,62,8,2,10,2,12,2,65,9,2,1,2,4,2,68,8,2,11,2,12,
        2,69,1,2,5,2,73,8,2,10,2,12,2,76,9,2,1,2,1,2,3,2,80,8,2,1,2,5,2,
        83,8,2,10,2,12,2,86,9,2,3,2,88,8,2,1,3,1,3,1,4,1,4,1,4,1,4,3,4,96,
        8,4,1,4,3,4,99,8,4,1,5,1,5,1,5,3,5,104,8,5,1,5,1,5,1,5,3,5,109,8,
        5,1,5,5,5,112,8,5,10,5,12,5,115,9,5,1,5,1,5,1,5,1,5,1,5,1,5,5,5,
        123,8,5,10,5,12,5,126,9,5,1,5,3,5,129,8,5,1,6,1,6,1,6,3,6,134,8,
        6,1,7,3,7,137,8,7,1,7,1,7,1,7,5,7,142,8,7,10,7,12,7,145,9,7,1,8,
        1,8,3,8,149,8,8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,
        3,8,163,8,8,1,8,3,8,166,8,8,1,8,5,8,169,8,8,10,8,12,8,172,9,8,1,
        8,1,8,1,8,1,8,1,8,1,8,1,8,3,8,181,8,8,1,8,1,8,1,8,3,8,186,8,8,1,
        8,1,8,1,8,1,8,1,8,3,8,193,8,8,1,9,1,9,1,9,1,9,5,9,199,8,9,10,9,12,
        9,202,9,9,3,9,204,8,9,1,10,1,10,1,11,1,11,1,12,1,12,1,13,1,13,1,
        13,0,0,14,0,2,4,6,8,10,12,14,16,18,20,22,24,26,0,6,2,0,3,3,13,13,
        2,0,27,28,34,34,2,0,14,14,36,36,1,0,4,7,2,0,29,29,32,32,1,0,1,2,
        237,0,28,1,0,0,0,2,54,1,0,0,0,4,87,1,0,0,0,6,89,1,0,0,0,8,91,1,0,
        0,0,10,128,1,0,0,0,12,133,1,0,0,0,14,136,1,0,0,0,16,192,1,0,0,0,
        18,203,1,0,0,0,20,205,1,0,0,0,22,207,1,0,0,0,24,209,1,0,0,0,26,211,
        1,0,0,0,28,29,3,26,13,0,29,32,3,2,1,0,30,31,5,26,0,0,31,33,3,2,1,
        0,32,30,1,0,0,0,32,33,1,0,0,0,33,34,1,0,0,0,34,35,5,37,0,0,35,1,
        1,0,0,0,36,37,3,8,4,0,37,38,3,4,2,0,38,55,1,0,0,0,39,40,5,21,0,0,
        40,41,3,2,1,0,41,42,5,22,0,0,42,43,3,2,1,0,43,55,1,0,0,0,44,45,3,
        8,4,0,45,50,3,4,2,0,46,47,5,21,0,0,47,48,3,2,1,0,48,49,5,22,0,0,
        49,51,1,0,0,0,50,46,1,0,0,0,50,51,1,0,0,0,51,52,1,0,0,0,52,53,3,
        2,1,0,53,55,1,0,0,0,54,36,1,0,0,0,54,39,1,0,0,0,54,44,1,0,0,0,55,
        3,1,0,0,0,56,57,5,11,0,0,57,59,5,19,0,0,58,56,1,0,0,0,58,59,1,0,
        0,0,59,63,1,0,0,0,60,62,3,16,8,0,61,60,1,0,0,0,62,65,1,0,0,0,63,
        61,1,0,0,0,63,64,1,0,0,0,64,67,1,0,0,0,65,63,1,0,0,0,66,68,3,6,3,
        0,67,66,1,0,0,0,68,69,1,0,0,0,69,67,1,0,0,0,69,70,1,0,0,0,70,74,
        1,0,0,0,71,73,3,16,8,0,72,71,1,0,0,0,73,76,1,0,0,0,74,72,1,0,0,0,
        74,75,1,0,0,0,75,79,1,0,0,0,76,74,1,0,0,0,77,78,5,18,0,0,78,80,5,
        15,0,0,79,77,1,0,0,0,79,80,1,0,0,0,80,88,1,0,0,0,81,83,3,16,8,0,
        82,81,1,0,0,0,83,86,1,0,0,0,84,82,1,0,0,0,84,85,1,0,0,0,85,88,1,
        0,0,0,86,84,1,0,0,0,87,58,1,0,0,0,87,84,1,0,0,0,88,5,1,0,0,0,89,
        90,7,0,0,0,90,7,1,0,0,0,91,95,5,20,0,0,92,93,3,18,9,0,93,94,3,20,
        10,0,94,96,1,0,0,0,95,92,1,0,0,0,95,96,1,0,0,0,96,98,1,0,0,0,97,
        99,3,18,9,0,98,97,1,0,0,0,98,99,1,0,0,0,99,9,1,0,0,0,100,101,5,33,
        0,0,101,103,5,25,0,0,102,104,3,24,12,0,103,102,1,0,0,0,103,104,1,
        0,0,0,104,105,1,0,0,0,105,113,5,16,0,0,106,108,5,18,0,0,107,109,
        3,24,12,0,108,107,1,0,0,0,108,109,1,0,0,0,109,110,1,0,0,0,110,112,
        5,16,0,0,111,106,1,0,0,0,112,115,1,0,0,0,113,111,1,0,0,0,113,114,
        1,0,0,0,114,116,1,0,0,0,115,113,1,0,0,0,116,129,5,26,0,0,117,118,
        5,29,0,0,118,119,5,25,0,0,119,124,5,16,0,0,120,121,5,18,0,0,121,
        123,5,16,0,0,122,120,1,0,0,0,123,126,1,0,0,0,124,122,1,0,0,0,124,
        125,1,0,0,0,125,127,1,0,0,0,126,124,1,0,0,0,127,129,5,26,0,0,128,
        100,1,0,0,0,128,117,1,0,0,0,129,11,1,0,0,0,130,134,5,13,0,0,131,
        134,3,22,11,0,132,134,5,8,0,0,133,130,1,0,0,0,133,131,1,0,0,0,133,
        132,1,0,0,0,134,13,1,0,0,0,135,137,7,1,0,0,136,135,1,0,0,0,136,137,
        1,0,0,0,137,138,1,0,0,0,138,139,5,4,0,0,139,143,5,16,0,0,140,142,
        3,10,5,0,141,140,1,0,0,0,142,145,1,0,0,0,143,141,1,0,0,0,143,144,
        1,0,0,0,144,15,1,0,0,0,145,143,1,0,0,0,146,147,5,16,0,0,147,149,
        5,18,0,0,148,146,1,0,0,0,148,149,1,0,0,0,149,150,1,0,0,0,150,151,
        5,16,0,0,151,152,5,19,0,0,152,153,5,9,0,0,153,193,5,19,0,0,154,155,
        5,16,0,0,155,156,5,19,0,0,156,157,3,22,11,0,157,158,5,19,0,0,158,
        159,3,12,6,0,159,160,5,19,0,0,160,193,1,0,0,0,161,163,5,19,0,0,162,
        161,1,0,0,0,162,163,1,0,0,0,163,165,1,0,0,0,164,166,5,16,0,0,165,
        164,1,0,0,0,165,166,1,0,0,0,166,170,1,0,0,0,167,169,3,22,11,0,168,
        167,1,0,0,0,169,172,1,0,0,0,170,168,1,0,0,0,170,171,1,0,0,0,171,
        173,1,0,0,0,172,170,1,0,0,0,173,193,3,12,6,0,174,175,5,36,0,0,175,
        176,5,16,0,0,176,193,5,8,0,0,177,185,5,16,0,0,178,179,5,18,0,0,179,
        181,5,16,0,0,180,178,1,0,0,0,180,181,1,0,0,0,181,182,1,0,0,0,182,
        186,5,30,0,0,183,186,5,31,0,0,184,186,3,14,7,0,185,180,1,0,0,0,185,
        183,1,0,0,0,185,184,1,0,0,0,186,193,1,0,0,0,187,193,5,10,0,0,188,
        189,5,11,0,0,189,193,5,19,0,0,190,191,5,19,0,0,191,193,5,12,0,0,
        192,148,1,0,0,0,192,154,1,0,0,0,192,162,1,0,0,0,192,174,1,0,0,0,
        192,177,1,0,0,0,192,187,1,0,0,0,192,188,1,0,0,0,192,190,1,0,0,0,
        193,17,1,0,0,0,194,204,5,36,0,0,195,200,5,16,0,0,196,197,5,35,0,
        0,197,199,5,16,0,0,198,196,1,0,0,0,199,202,1,0,0,0,200,198,1,0,0,
        0,200,201,1,0,0,0,201,204,1,0,0,0,202,200,1,0,0,0,203,194,1,0,0,
        0,203,195,1,0,0,0,204,19,1,0,0,0,205,206,7,2,0,0,206,21,1,0,0,0,
        207,208,7,3,0,0,208,23,1,0,0,0,209,210,7,4,0,0,210,25,1,0,0,0,211,
        212,7,5,0,0,212,27,1,0,0,0,29,32,50,54,58,63,69,74,79,84,87,95,98,
        103,108,113,124,128,133,136,143,148,162,165,170,180,185,192,200,
        203
    ]

class GWBParser ( Parser ):

    grammarFileName = "GWB.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'freeEnd'", "'redEnd'", "<INVALID>", 
                     "'C'", "'N'", "'O'", "'P'", "<INVALID>", "'Anhydro'", 
                     "'0d'", "<INVALID>", "<INVALID>", "<INVALID>", "<INVALID>", 
                     "<INVALID>", "<INVALID>", "';'", "','", "'-'", "'--'", 
                     "'('", "')'", "'['", "']'", "'{'", "'}'", "'a'", "'ai'", 
                     "'c'", "'d'", "'e'", "'t'", "'='", "'i'", "'/'", "'?'", 
                     "'$'" ]

    symbolicNames = [ "<INVALID>", "FREEEND", "REDEND", "SAC", "CARBON", 
                      "NITROGEN", "OXYGEN", "PHOSPHOR", "FG", "ANHYDRO", 
                      "HEAD", "HEADD", "END", "COUNT", "TYPE", "RING", "NUM", 
                      "SEMICOLON", "COLON", "DASH", "DOUBLEDASH", "LPAR", 
                      "RPAR", "LBRACE", "RBRACE", "LBRACK", "RBRACK", "A", 
                      "AI", "C", "D", "E", "T", "EQ", "I", "SLASH", "QMARK", 
                      "DOLLAR" ]

    RULE_start = 0
    RULE_branch = 1
    RULE_deriv = 2
    RULE_saci = 3
    RULE_con = 4
    RULE_add = 5
    RULE_fgi = 6
    RULE_carb = 7
    RULE_modi = 8
    RULE_qnum = 9
    RULE_typi = 10
    RULE_bridge = 11
    RULE_ct = 12
    RULE_bos = 13

    ruleNames =  [ "start", "branch", "deriv", "saci", "con", "add", "fgi", 
                   "carb", "modi", "qnum", "typi", "bridge", "ct", "bos" ]

    EOF = Token.EOF
    FREEEND=1
    REDEND=2
    SAC=3
    CARBON=4
    NITROGEN=5
    OXYGEN=6
    PHOSPHOR=7
    FG=8
    ANHYDRO=9
    HEAD=10
    HEADD=11
    END=12
    COUNT=13
    TYPE=14
    RING=15
    NUM=16
    SEMICOLON=17
    COLON=18
    DASH=19
    DOUBLEDASH=20
    LPAR=21
    RPAR=22
    LBRACE=23
    RBRACE=24
    LBRACK=25
    RBRACK=26
    A=27
    AI=28
    C=29
    D=30
    E=31
    T=32
    EQ=33
    I=34
    SLASH=35
    QMARK=36
    DOLLAR=37

    def __init__(self, input:TokenStream, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.13.2")
        self._interp = ParserATNSimulator(self, self.atn, self.decisionsToDFA, self.sharedContextCache)
        self._predicates = None




    class StartContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def bos(self):
            return self.getTypedRuleContext(GWBParser.BosContext,0)


        def branch(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GWBParser.BranchContext)
            else:
                return self.getTypedRuleContext(GWBParser.BranchContext,i)


        def DOLLAR(self):
            return self.getToken(GWBParser.DOLLAR, 0)

        def RBRACK(self):
            return self.getToken(GWBParser.RBRACK, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_start

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterStart" ):
                listener.enterStart(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitStart" ):
                listener.exitStart(self)




    def start(self):

        localctx = GWBParser.StartContext(self, self._ctx, self.state)
        self.enterRule(localctx, 0, self.RULE_start)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 28
            self.bos()
            self.state = 29
            self.branch()
            self.state = 32
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            if _la==26:
                self.state = 30
                self.match(GWBParser.RBRACK)
                self.state = 31
                self.branch()


            self.state = 34
            self.match(GWBParser.DOLLAR)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class BranchContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def con(self):
            return self.getTypedRuleContext(GWBParser.ConContext,0)


        def deriv(self):
            return self.getTypedRuleContext(GWBParser.DerivContext,0)


        def LPAR(self):
            return self.getToken(GWBParser.LPAR, 0)

        def branch(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GWBParser.BranchContext)
            else:
                return self.getTypedRuleContext(GWBParser.BranchContext,i)


        def RPAR(self):
            return self.getToken(GWBParser.RPAR, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_branch

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterBranch" ):
                listener.enterBranch(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitBranch" ):
                listener.exitBranch(self)




    def branch(self):

        localctx = GWBParser.BranchContext(self, self._ctx, self.state)
        self.enterRule(localctx, 2, self.RULE_branch)
        try:
            self.state = 54
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,2,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 36
                self.con()
                self.state = 37
                self.deriv()
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 39
                self.match(GWBParser.LPAR)
                self.state = 40
                self.branch()
                self.state = 41
                self.match(GWBParser.RPAR)
                self.state = 42
                self.branch()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 44
                self.con()
                self.state = 45
                self.deriv()
                self.state = 50
                self._errHandler.sync(self)
                la_ = self._interp.adaptivePredict(self._input,1,self._ctx)
                if la_ == 1:
                    self.state = 46
                    self.match(GWBParser.LPAR)
                    self.state = 47
                    self.branch()
                    self.state = 48
                    self.match(GWBParser.RPAR)


                self.state = 52
                self.branch()
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class DerivContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def HEADD(self):
            return self.getToken(GWBParser.HEADD, 0)

        def DASH(self):
            return self.getToken(GWBParser.DASH, 0)

        def modi(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GWBParser.ModiContext)
            else:
                return self.getTypedRuleContext(GWBParser.ModiContext,i)


        def saci(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GWBParser.SaciContext)
            else:
                return self.getTypedRuleContext(GWBParser.SaciContext,i)


        def COLON(self):
            return self.getToken(GWBParser.COLON, 0)

        def RING(self):
            return self.getToken(GWBParser.RING, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_deriv

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterDeriv" ):
                listener.enterDeriv(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitDeriv" ):
                listener.exitDeriv(self)




    def deriv(self):

        localctx = GWBParser.DerivContext(self, self._ctx, self.state)
        self.enterRule(localctx, 4, self.RULE_deriv)
        self._la = 0 # Token type
        try:
            self.state = 87
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,9,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 58
                self._errHandler.sync(self)
                la_ = self._interp.adaptivePredict(self._input,3,self._ctx)
                if la_ == 1:
                    self.state = 56
                    self.match(GWBParser.HEADD)
                    self.state = 57
                    self.match(GWBParser.DASH)


                self.state = 63
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,4,self._ctx)
                while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                    if _alt==1:
                        self.state = 60
                        self.modi() 
                    self.state = 65
                    self._errHandler.sync(self)
                    _alt = self._interp.adaptivePredict(self._input,4,self._ctx)

                self.state = 67 
                self._errHandler.sync(self)
                _alt = 1
                while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                    if _alt == 1:
                        self.state = 66
                        self.saci()

                    else:
                        raise NoViableAltException(self)
                    self.state = 69 
                    self._errHandler.sync(self)
                    _alt = self._interp.adaptivePredict(self._input,5,self._ctx)

                self.state = 74
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while (((_la) & ~0x3f) == 0 and ((1 << _la) & 68720078320) != 0):
                    self.state = 71
                    self.modi()
                    self.state = 76
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                self.state = 79
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==18:
                    self.state = 77
                    self.match(GWBParser.COLON)
                    self.state = 78
                    self.match(GWBParser.RING)


                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 84
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while (((_la) & ~0x3f) == 0 and ((1 << _la) & 68720078320) != 0):
                    self.state = 81
                    self.modi()
                    self.state = 86
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class SaciContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def COUNT(self):
            return self.getToken(GWBParser.COUNT, 0)

        def SAC(self):
            return self.getToken(GWBParser.SAC, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_saci

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSaci" ):
                listener.enterSaci(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSaci" ):
                listener.exitSaci(self)




    def saci(self):

        localctx = GWBParser.SaciContext(self, self._ctx, self.state)
        self.enterRule(localctx, 6, self.RULE_saci)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 89
            _la = self._input.LA(1)
            if not(_la==3 or _la==13):
                self._errHandler.recoverInline(self)
            else:
                self._errHandler.reportMatch(self)
                self.consume()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class ConContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def DOUBLEDASH(self):
            return self.getToken(GWBParser.DOUBLEDASH, 0)

        def qnum(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GWBParser.QnumContext)
            else:
                return self.getTypedRuleContext(GWBParser.QnumContext,i)


        def typi(self):
            return self.getTypedRuleContext(GWBParser.TypiContext,0)


        def getRuleIndex(self):
            return GWBParser.RULE_con

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterCon" ):
                listener.enterCon(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitCon" ):
                listener.exitCon(self)




    def con(self):

        localctx = GWBParser.ConContext(self, self._ctx, self.state)
        self.enterRule(localctx, 8, self.RULE_con)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 91
            self.match(GWBParser.DOUBLEDASH)
            self.state = 95
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,10,self._ctx)
            if la_ == 1:
                self.state = 92
                self.qnum()
                self.state = 93
                self.typi()


            self.state = 98
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,11,self._ctx)
            if la_ == 1:
                self.state = 97
                self.qnum()


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class AddContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def EQ(self):
            return self.getToken(GWBParser.EQ, 0)

        def LBRACK(self):
            return self.getToken(GWBParser.LBRACK, 0)

        def NUM(self, i:int=None):
            if i is None:
                return self.getTokens(GWBParser.NUM)
            else:
                return self.getToken(GWBParser.NUM, i)

        def RBRACK(self):
            return self.getToken(GWBParser.RBRACK, 0)

        def ct(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GWBParser.CtContext)
            else:
                return self.getTypedRuleContext(GWBParser.CtContext,i)


        def COLON(self, i:int=None):
            if i is None:
                return self.getTokens(GWBParser.COLON)
            else:
                return self.getToken(GWBParser.COLON, i)

        def C(self):
            return self.getToken(GWBParser.C, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_add

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterAdd" ):
                listener.enterAdd(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitAdd" ):
                listener.exitAdd(self)




    def add(self):

        localctx = GWBParser.AddContext(self, self._ctx, self.state)
        self.enterRule(localctx, 10, self.RULE_add)
        self._la = 0 # Token type
        try:
            self.state = 128
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [33]:
                self.enterOuterAlt(localctx, 1)
                self.state = 100
                self.match(GWBParser.EQ)
                self.state = 101
                self.match(GWBParser.LBRACK)
                self.state = 103
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==29 or _la==32:
                    self.state = 102
                    self.ct()


                self.state = 105
                self.match(GWBParser.NUM)
                self.state = 113
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==18:
                    self.state = 106
                    self.match(GWBParser.COLON)
                    self.state = 108
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)
                    if _la==29 or _la==32:
                        self.state = 107
                        self.ct()


                    self.state = 110
                    self.match(GWBParser.NUM)
                    self.state = 115
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                self.state = 116
                self.match(GWBParser.RBRACK)
                pass
            elif token in [29]:
                self.enterOuterAlt(localctx, 2)
                self.state = 117
                self.match(GWBParser.C)
                self.state = 118
                self.match(GWBParser.LBRACK)
                self.state = 119
                self.match(GWBParser.NUM)
                self.state = 124
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==18:
                    self.state = 120
                    self.match(GWBParser.COLON)
                    self.state = 121
                    self.match(GWBParser.NUM)
                    self.state = 126
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                self.state = 127
                self.match(GWBParser.RBRACK)
                pass
            else:
                raise NoViableAltException(self)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class FgiContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def COUNT(self):
            return self.getToken(GWBParser.COUNT, 0)

        def bridge(self):
            return self.getTypedRuleContext(GWBParser.BridgeContext,0)


        def FG(self):
            return self.getToken(GWBParser.FG, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_fgi

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterFgi" ):
                listener.enterFgi(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitFgi" ):
                listener.exitFgi(self)




    def fgi(self):

        localctx = GWBParser.FgiContext(self, self._ctx, self.state)
        self.enterRule(localctx, 12, self.RULE_fgi)
        try:
            self.state = 133
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [13]:
                self.enterOuterAlt(localctx, 1)
                self.state = 130
                self.match(GWBParser.COUNT)
                pass
            elif token in [4, 5, 6, 7]:
                self.enterOuterAlt(localctx, 2)
                self.state = 131
                self.bridge()
                pass
            elif token in [8]:
                self.enterOuterAlt(localctx, 3)
                self.state = 132
                self.match(GWBParser.FG)
                pass
            else:
                raise NoViableAltException(self)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class CarbContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def CARBON(self):
            return self.getToken(GWBParser.CARBON, 0)

        def NUM(self):
            return self.getToken(GWBParser.NUM, 0)

        def add(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GWBParser.AddContext)
            else:
                return self.getTypedRuleContext(GWBParser.AddContext,i)


        def AI(self):
            return self.getToken(GWBParser.AI, 0)

        def A(self):
            return self.getToken(GWBParser.A, 0)

        def I(self):
            return self.getToken(GWBParser.I, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_carb

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterCarb" ):
                listener.enterCarb(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitCarb" ):
                listener.exitCarb(self)




    def carb(self):

        localctx = GWBParser.CarbContext(self, self._ctx, self.state)
        self.enterRule(localctx, 14, self.RULE_carb)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 136
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            if (((_la) & ~0x3f) == 0 and ((1 << _la) & 17582522368) != 0):
                self.state = 135
                _la = self._input.LA(1)
                if not((((_la) & ~0x3f) == 0 and ((1 << _la) & 17582522368) != 0)):
                    self._errHandler.recoverInline(self)
                else:
                    self._errHandler.reportMatch(self)
                    self.consume()


            self.state = 138
            self.match(GWBParser.CARBON)
            self.state = 139
            self.match(GWBParser.NUM)
            self.state = 143
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==29 or _la==33:
                self.state = 140
                self.add()
                self.state = 145
                self._errHandler.sync(self)
                _la = self._input.LA(1)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class ModiContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def NUM(self, i:int=None):
            if i is None:
                return self.getTokens(GWBParser.NUM)
            else:
                return self.getToken(GWBParser.NUM, i)

        def DASH(self, i:int=None):
            if i is None:
                return self.getTokens(GWBParser.DASH)
            else:
                return self.getToken(GWBParser.DASH, i)

        def ANHYDRO(self):
            return self.getToken(GWBParser.ANHYDRO, 0)

        def COLON(self):
            return self.getToken(GWBParser.COLON, 0)

        def bridge(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(GWBParser.BridgeContext)
            else:
                return self.getTypedRuleContext(GWBParser.BridgeContext,i)


        def fgi(self):
            return self.getTypedRuleContext(GWBParser.FgiContext,0)


        def QMARK(self):
            return self.getToken(GWBParser.QMARK, 0)

        def FG(self):
            return self.getToken(GWBParser.FG, 0)

        def D(self):
            return self.getToken(GWBParser.D, 0)

        def E(self):
            return self.getToken(GWBParser.E, 0)

        def carb(self):
            return self.getTypedRuleContext(GWBParser.CarbContext,0)


        def HEAD(self):
            return self.getToken(GWBParser.HEAD, 0)

        def HEADD(self):
            return self.getToken(GWBParser.HEADD, 0)

        def END(self):
            return self.getToken(GWBParser.END, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_modi

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterModi" ):
                listener.enterModi(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitModi" ):
                listener.exitModi(self)




    def modi(self):

        localctx = GWBParser.ModiContext(self, self._ctx, self.state)
        self.enterRule(localctx, 16, self.RULE_modi)
        self._la = 0 # Token type
        try:
            self.state = 192
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,26,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 148
                self._errHandler.sync(self)
                la_ = self._interp.adaptivePredict(self._input,20,self._ctx)
                if la_ == 1:
                    self.state = 146
                    self.match(GWBParser.NUM)
                    self.state = 147
                    self.match(GWBParser.COLON)


                self.state = 150
                self.match(GWBParser.NUM)
                self.state = 151
                self.match(GWBParser.DASH)
                self.state = 152
                self.match(GWBParser.ANHYDRO)
                self.state = 153
                self.match(GWBParser.DASH)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 154
                self.match(GWBParser.NUM)
                self.state = 155
                self.match(GWBParser.DASH)
                self.state = 156
                self.bridge()
                self.state = 157
                self.match(GWBParser.DASH)
                self.state = 158
                self.fgi()
                self.state = 159
                self.match(GWBParser.DASH)
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 162
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==19:
                    self.state = 161
                    self.match(GWBParser.DASH)


                self.state = 165
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==16:
                    self.state = 164
                    self.match(GWBParser.NUM)


                self.state = 170
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,23,self._ctx)
                while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                    if _alt==1:
                        self.state = 167
                        self.bridge() 
                    self.state = 172
                    self._errHandler.sync(self)
                    _alt = self._interp.adaptivePredict(self._input,23,self._ctx)

                self.state = 173
                self.fgi()
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 174
                self.match(GWBParser.QMARK)
                self.state = 175
                self.match(GWBParser.NUM)
                self.state = 176
                self.match(GWBParser.FG)
                pass

            elif la_ == 5:
                self.enterOuterAlt(localctx, 5)
                self.state = 177
                self.match(GWBParser.NUM)
                self.state = 185
                self._errHandler.sync(self)
                token = self._input.LA(1)
                if token in [18, 30]:
                    self.state = 180
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)
                    if _la==18:
                        self.state = 178
                        self.match(GWBParser.COLON)
                        self.state = 179
                        self.match(GWBParser.NUM)


                    self.state = 182
                    self.match(GWBParser.D)
                    pass
                elif token in [31]:
                    self.state = 183
                    self.match(GWBParser.E)
                    pass
                elif token in [4, 27, 28, 34]:
                    self.state = 184
                    self.carb()
                    pass
                else:
                    raise NoViableAltException(self)

                pass

            elif la_ == 6:
                self.enterOuterAlt(localctx, 6)
                self.state = 187
                self.match(GWBParser.HEAD)
                pass

            elif la_ == 7:
                self.enterOuterAlt(localctx, 7)
                self.state = 188
                self.match(GWBParser.HEADD)
                self.state = 189
                self.match(GWBParser.DASH)
                pass

            elif la_ == 8:
                self.enterOuterAlt(localctx, 8)
                self.state = 190
                self.match(GWBParser.DASH)
                self.state = 191
                self.match(GWBParser.END)
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class QnumContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def QMARK(self):
            return self.getToken(GWBParser.QMARK, 0)

        def NUM(self, i:int=None):
            if i is None:
                return self.getTokens(GWBParser.NUM)
            else:
                return self.getToken(GWBParser.NUM, i)

        def SLASH(self, i:int=None):
            if i is None:
                return self.getTokens(GWBParser.SLASH)
            else:
                return self.getToken(GWBParser.SLASH, i)

        def getRuleIndex(self):
            return GWBParser.RULE_qnum

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterQnum" ):
                listener.enterQnum(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitQnum" ):
                listener.exitQnum(self)




    def qnum(self):

        localctx = GWBParser.QnumContext(self, self._ctx, self.state)
        self.enterRule(localctx, 18, self.RULE_qnum)
        self._la = 0 # Token type
        try:
            self.state = 203
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [36]:
                self.enterOuterAlt(localctx, 1)
                self.state = 194
                self.match(GWBParser.QMARK)
                pass
            elif token in [16]:
                self.enterOuterAlt(localctx, 2)
                self.state = 195
                self.match(GWBParser.NUM)
                self.state = 200
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==35:
                    self.state = 196
                    self.match(GWBParser.SLASH)
                    self.state = 197
                    self.match(GWBParser.NUM)
                    self.state = 202
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                pass
            else:
                raise NoViableAltException(self)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class TypiContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def TYPE(self):
            return self.getToken(GWBParser.TYPE, 0)

        def QMARK(self):
            return self.getToken(GWBParser.QMARK, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_typi

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterTypi" ):
                listener.enterTypi(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitTypi" ):
                listener.exitTypi(self)




    def typi(self):

        localctx = GWBParser.TypiContext(self, self._ctx, self.state)
        self.enterRule(localctx, 20, self.RULE_typi)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 205
            _la = self._input.LA(1)
            if not(_la==14 or _la==36):
                self._errHandler.recoverInline(self)
            else:
                self._errHandler.reportMatch(self)
                self.consume()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class BridgeContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def CARBON(self):
            return self.getToken(GWBParser.CARBON, 0)

        def NITROGEN(self):
            return self.getToken(GWBParser.NITROGEN, 0)

        def OXYGEN(self):
            return self.getToken(GWBParser.OXYGEN, 0)

        def PHOSPHOR(self):
            return self.getToken(GWBParser.PHOSPHOR, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_bridge

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterBridge" ):
                listener.enterBridge(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitBridge" ):
                listener.exitBridge(self)




    def bridge(self):

        localctx = GWBParser.BridgeContext(self, self._ctx, self.state)
        self.enterRule(localctx, 22, self.RULE_bridge)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 207
            _la = self._input.LA(1)
            if not((((_la) & ~0x3f) == 0 and ((1 << _la) & 240) != 0)):
                self._errHandler.recoverInline(self)
            else:
                self._errHandler.reportMatch(self)
                self.consume()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class CtContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def C(self):
            return self.getToken(GWBParser.C, 0)

        def T(self):
            return self.getToken(GWBParser.T, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_ct

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterCt" ):
                listener.enterCt(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitCt" ):
                listener.exitCt(self)




    def ct(self):

        localctx = GWBParser.CtContext(self, self._ctx, self.state)
        self.enterRule(localctx, 24, self.RULE_ct)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 209
            _la = self._input.LA(1)
            if not(_la==29 or _la==32):
                self._errHandler.recoverInline(self)
            else:
                self._errHandler.reportMatch(self)
                self.consume()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class BosContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def FREEEND(self):
            return self.getToken(GWBParser.FREEEND, 0)

        def REDEND(self):
            return self.getToken(GWBParser.REDEND, 0)

        def getRuleIndex(self):
            return GWBParser.RULE_bos

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterBos" ):
                listener.enterBos(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitBos" ):
                listener.exitBos(self)




    def bos(self):

        localctx = GWBParser.BosContext(self, self._ctx, self.state)
        self.enterRule(localctx, 26, self.RULE_bos)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 211
            _la = self._input.LA(1)
            if not(_la==1 or _la==2):
                self._errHandler.recoverInline(self)
            else:
                self._errHandler.reportMatch(self)
                self.consume()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





