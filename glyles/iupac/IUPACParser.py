# Generated from IUPAC.g4 by ANTLR 4.13.2
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
        4,1,41,263,2,0,7,0,2,1,7,1,2,2,7,2,2,3,7,3,2,4,7,4,2,5,7,5,2,6,7,
        6,2,7,7,7,2,8,7,8,2,9,7,9,2,10,7,10,2,11,7,11,2,12,7,12,2,13,7,13,
        1,0,1,0,1,0,1,0,1,0,5,0,34,8,0,10,0,12,0,37,9,0,1,0,1,0,1,0,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,55,8,1,1,2,1,
        2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,
        2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,
        2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,3,2,98,8,2,1,3,5,3,101,8,3,10,3,12,
        3,104,9,3,1,3,4,3,107,8,3,11,3,12,3,108,1,3,5,3,112,8,3,10,3,12,
        3,115,9,3,1,3,3,3,118,8,3,1,3,5,3,121,8,3,10,3,12,3,124,9,3,1,3,
        3,3,127,8,3,1,4,1,4,1,5,1,5,1,5,1,5,1,5,1,5,1,5,1,5,1,5,1,5,1,5,
        1,5,1,5,1,5,1,5,1,5,1,5,1,5,1,5,1,5,1,5,3,5,152,8,5,1,6,1,6,1,6,
        3,6,157,8,6,1,6,1,6,1,6,3,6,162,8,6,1,6,5,6,165,8,6,10,6,12,6,168,
        9,6,1,6,1,6,1,6,1,6,1,6,1,6,5,6,176,8,6,10,6,12,6,179,9,6,1,6,3,
        6,182,8,6,1,7,1,7,1,7,3,7,187,8,7,1,8,3,8,190,8,8,1,8,1,8,1,8,5,
        8,195,8,8,10,8,12,8,198,9,8,1,9,1,9,3,9,202,8,9,1,9,1,9,1,9,1,9,
        1,9,1,9,1,9,1,9,1,9,1,9,1,9,1,9,3,9,216,8,9,1,9,3,9,219,8,9,1,9,
        5,9,222,8,9,10,9,12,9,225,9,9,1,9,1,9,1,9,1,9,3,9,231,8,9,1,9,1,
        9,1,9,1,9,3,9,237,8,9,1,9,1,9,1,9,1,9,1,9,3,9,244,8,9,1,10,1,10,
        1,10,1,10,5,10,250,8,10,10,10,12,10,253,9,10,3,10,255,8,10,1,11,
        1,11,1,12,1,12,1,13,1,13,1,13,0,0,14,0,2,4,6,8,10,12,14,16,18,20,
        22,24,26,0,5,2,0,3,3,13,13,2,0,28,29,35,35,2,0,14,14,41,41,1,0,4,
        7,2,0,30,30,33,33,291,0,28,1,0,0,0,2,54,1,0,0,0,4,97,1,0,0,0,6,102,
        1,0,0,0,8,128,1,0,0,0,10,151,1,0,0,0,12,181,1,0,0,0,14,186,1,0,0,
        0,16,189,1,0,0,0,18,243,1,0,0,0,20,254,1,0,0,0,22,256,1,0,0,0,24,
        258,1,0,0,0,26,260,1,0,0,0,28,35,5,39,0,0,29,30,5,26,0,0,30,31,3,
        4,2,0,31,32,5,27,0,0,32,34,1,0,0,0,33,29,1,0,0,0,34,37,1,0,0,0,35,
        33,1,0,0,0,35,36,1,0,0,0,36,38,1,0,0,0,37,35,1,0,0,0,38,39,3,2,1,
        0,39,40,5,39,0,0,40,1,1,0,0,0,41,42,3,4,2,0,42,43,3,6,3,0,43,44,
        5,40,0,0,44,45,5,14,0,0,45,55,1,0,0,0,46,47,3,4,2,0,47,48,3,6,3,
        0,48,55,1,0,0,0,49,50,3,6,3,0,50,51,5,40,0,0,51,52,5,14,0,0,52,55,
        1,0,0,0,53,55,3,6,3,0,54,41,1,0,0,0,54,46,1,0,0,0,54,49,1,0,0,0,
        54,53,1,0,0,0,55,3,1,0,0,0,56,57,3,6,3,0,57,58,3,10,5,0,58,98,1,
        0,0,0,59,60,3,6,3,0,60,61,3,10,5,0,61,62,3,4,2,0,62,98,1,0,0,0,63,
        64,5,24,0,0,64,65,3,4,2,0,65,66,5,25,0,0,66,98,1,0,0,0,67,68,3,6,
        3,0,68,69,3,10,5,0,69,70,5,24,0,0,70,71,3,4,2,0,71,72,5,25,0,0,72,
        73,3,4,2,0,73,98,1,0,0,0,74,75,3,6,3,0,75,76,3,10,5,0,76,77,5,24,
        0,0,77,78,3,4,2,0,78,79,5,25,0,0,79,80,5,24,0,0,80,81,3,4,2,0,81,
        82,5,25,0,0,82,83,3,4,2,0,83,98,1,0,0,0,84,85,3,6,3,0,85,86,3,10,
        5,0,86,87,5,24,0,0,87,88,3,4,2,0,88,89,5,25,0,0,89,90,5,24,0,0,90,
        91,3,4,2,0,91,92,5,25,0,0,92,93,5,24,0,0,93,94,3,4,2,0,94,95,5,25,
        0,0,95,96,3,4,2,0,96,98,1,0,0,0,97,56,1,0,0,0,97,59,1,0,0,0,97,63,
        1,0,0,0,97,67,1,0,0,0,97,74,1,0,0,0,97,84,1,0,0,0,98,5,1,0,0,0,99,
        101,3,18,9,0,100,99,1,0,0,0,101,104,1,0,0,0,102,100,1,0,0,0,102,
        103,1,0,0,0,103,106,1,0,0,0,104,102,1,0,0,0,105,107,3,8,4,0,106,
        105,1,0,0,0,107,108,1,0,0,0,108,106,1,0,0,0,108,109,1,0,0,0,109,
        113,1,0,0,0,110,112,3,18,9,0,111,110,1,0,0,0,112,115,1,0,0,0,113,
        111,1,0,0,0,113,114,1,0,0,0,114,117,1,0,0,0,115,113,1,0,0,0,116,
        118,5,15,0,0,117,116,1,0,0,0,117,118,1,0,0,0,118,122,1,0,0,0,119,
        121,3,18,9,0,120,119,1,0,0,0,121,124,1,0,0,0,122,120,1,0,0,0,122,
        123,1,0,0,0,123,126,1,0,0,0,124,122,1,0,0,0,125,127,5,14,0,0,126,
        125,1,0,0,0,126,127,1,0,0,0,127,7,1,0,0,0,128,129,7,0,0,0,129,9,
        1,0,0,0,130,131,5,22,0,0,131,132,3,22,11,0,132,133,5,16,0,0,133,
        134,5,21,0,0,134,135,3,20,10,0,135,136,5,23,0,0,136,152,1,0,0,0,
        137,138,5,22,0,0,138,139,5,16,0,0,139,140,5,21,0,0,140,141,3,20,
        10,0,141,142,5,23,0,0,142,152,1,0,0,0,143,144,3,22,11,0,144,145,
        5,16,0,0,145,146,5,21,0,0,146,147,3,20,10,0,147,152,1,0,0,0,148,
        149,3,22,11,0,149,150,5,16,0,0,150,152,1,0,0,0,151,130,1,0,0,0,151,
        137,1,0,0,0,151,143,1,0,0,0,151,148,1,0,0,0,152,11,1,0,0,0,153,154,
        5,34,0,0,154,156,5,26,0,0,155,157,3,26,13,0,156,155,1,0,0,0,156,
        157,1,0,0,0,157,158,1,0,0,0,158,166,5,16,0,0,159,161,5,18,0,0,160,
        162,3,26,13,0,161,160,1,0,0,0,161,162,1,0,0,0,162,163,1,0,0,0,163,
        165,5,16,0,0,164,159,1,0,0,0,165,168,1,0,0,0,166,164,1,0,0,0,166,
        167,1,0,0,0,167,169,1,0,0,0,168,166,1,0,0,0,169,182,5,27,0,0,170,
        171,5,30,0,0,171,172,5,26,0,0,172,177,5,16,0,0,173,174,5,18,0,0,
        174,176,5,16,0,0,175,173,1,0,0,0,176,179,1,0,0,0,177,175,1,0,0,0,
        177,178,1,0,0,0,178,180,1,0,0,0,179,177,1,0,0,0,180,182,5,27,0,0,
        181,153,1,0,0,0,181,170,1,0,0,0,182,13,1,0,0,0,183,187,5,13,0,0,
        184,187,3,24,12,0,185,187,5,8,0,0,186,183,1,0,0,0,186,184,1,0,0,
        0,186,185,1,0,0,0,187,15,1,0,0,0,188,190,7,1,0,0,189,188,1,0,0,0,
        189,190,1,0,0,0,190,191,1,0,0,0,191,192,5,4,0,0,192,196,5,16,0,0,
        193,195,3,12,6,0,194,193,1,0,0,0,195,198,1,0,0,0,196,194,1,0,0,0,
        196,197,1,0,0,0,197,17,1,0,0,0,198,196,1,0,0,0,199,200,5,16,0,0,
        200,202,5,18,0,0,201,199,1,0,0,0,201,202,1,0,0,0,202,203,1,0,0,0,
        203,204,5,16,0,0,204,205,5,21,0,0,205,206,5,9,0,0,206,244,5,21,0,
        0,207,208,5,16,0,0,208,209,5,21,0,0,209,210,3,24,12,0,210,211,5,
        21,0,0,211,212,3,14,7,0,212,213,5,21,0,0,213,244,1,0,0,0,214,216,
        5,21,0,0,215,214,1,0,0,0,215,216,1,0,0,0,216,218,1,0,0,0,217,219,
        5,16,0,0,218,217,1,0,0,0,218,219,1,0,0,0,219,223,1,0,0,0,220,222,
        3,24,12,0,221,220,1,0,0,0,222,225,1,0,0,0,223,221,1,0,0,0,223,224,
        1,0,0,0,224,226,1,0,0,0,225,223,1,0,0,0,226,244,3,14,7,0,227,236,
        5,16,0,0,228,229,5,18,0,0,229,231,5,16,0,0,230,228,1,0,0,0,230,231,
        1,0,0,0,231,232,1,0,0,0,232,237,5,31,0,0,233,237,5,32,0,0,234,237,
        5,36,0,0,235,237,3,16,8,0,236,230,1,0,0,0,236,233,1,0,0,0,236,234,
        1,0,0,0,236,235,1,0,0,0,237,244,1,0,0,0,238,244,5,10,0,0,239,240,
        5,11,0,0,240,244,5,21,0,0,241,242,5,21,0,0,242,244,5,12,0,0,243,
        201,1,0,0,0,243,207,1,0,0,0,243,215,1,0,0,0,243,227,1,0,0,0,243,
        238,1,0,0,0,243,239,1,0,0,0,243,241,1,0,0,0,244,19,1,0,0,0,245,255,
        5,41,0,0,246,251,5,16,0,0,247,248,5,37,0,0,248,250,5,16,0,0,249,
        247,1,0,0,0,250,253,1,0,0,0,251,249,1,0,0,0,251,252,1,0,0,0,252,
        255,1,0,0,0,253,251,1,0,0,0,254,245,1,0,0,0,254,246,1,0,0,0,255,
        21,1,0,0,0,256,257,7,2,0,0,257,23,1,0,0,0,258,259,7,3,0,0,259,25,
        1,0,0,0,260,261,7,4,0,0,261,27,1,0,0,0,27,35,54,97,102,108,113,117,
        122,126,151,156,161,166,177,181,186,189,196,201,215,218,223,230,
        236,243,251,254
    ]

class IUPACParser ( Parser ):

    grammarFileName = "IUPAC.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'freeEnd'", "'redEnd'", "<INVALID>", 
                     "'C'", "'N'", "'O'", "'P'", "<INVALID>", "'Anhydro'", 
                     "'0d'", "<INVALID>", "<INVALID>", "<INVALID>", "<INVALID>", 
                     "<INVALID>", "<INVALID>", "'@'", "','", "';'", "'--'", 
                     "'-'", "'('", "')'", "'['", "']'", "'{'", "'}'", "'a'", 
                     "'ai'", "'c'", "'d'", "'e'", "'t'", "'='", "'i'", "'u'", 
                     "'/'", "'$'", "'#'", "' '", "'?'" ]

    symbolicNames = [ "<INVALID>", "FREEEND", "REDEND", "SAC", "CARBON", 
                      "NITROGEN", "OXYGEN", "PHOSPHOR", "FG", "ANHYDRO", 
                      "HEAD", "HEADD", "END", "COUNT", "TYPE", "RING", "NUM", 
                      "AT", "COMMA", "SEMICOLON", "DOUBLEDASH", "DASH", 
                      "LPAR", "RPAR", "LBRACE", "RBRACE", "LBRACK", "RBRACK", 
                      "A", "AI", "C", "D", "E", "T", "EQ", "I", "U", "SLASH", 
                      "DOLLAR", "HASH", "SPACE", "QMARK" ]

    RULE_start = 0
    RULE_begin = 1
    RULE_branch = 2
    RULE_deriv = 3
    RULE_saci = 4
    RULE_con = 5
    RULE_add = 6
    RULE_fgi = 7
    RULE_carb = 8
    RULE_modi = 9
    RULE_qnum = 10
    RULE_typi = 11
    RULE_bridge = 12
    RULE_ct = 13

    ruleNames =  [ "start", "begin", "branch", "deriv", "saci", "con", "add", 
                   "fgi", "carb", "modi", "qnum", "typi", "bridge", "ct" ]

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
    AT=17
    COMMA=18
    SEMICOLON=19
    DOUBLEDASH=20
    DASH=21
    LPAR=22
    RPAR=23
    LBRACE=24
    RBRACE=25
    LBRACK=26
    RBRACK=27
    A=28
    AI=29
    C=30
    D=31
    E=32
    T=33
    EQ=34
    I=35
    U=36
    SLASH=37
    DOLLAR=38
    HASH=39
    SPACE=40
    QMARK=41

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

        def HASH(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.HASH)
            else:
                return self.getToken(IUPACParser.HASH, i)

        def begin(self):
            return self.getTypedRuleContext(IUPACParser.BeginContext,0)


        def LBRACK(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.LBRACK)
            else:
                return self.getToken(IUPACParser.LBRACK, i)

        def branch(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(IUPACParser.BranchContext)
            else:
                return self.getTypedRuleContext(IUPACParser.BranchContext,i)


        def RBRACK(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.RBRACK)
            else:
                return self.getToken(IUPACParser.RBRACK, i)

        def getRuleIndex(self):
            return IUPACParser.RULE_start

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterStart" ):
                listener.enterStart(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitStart" ):
                listener.exitStart(self)




    def start(self):

        localctx = IUPACParser.StartContext(self, self._ctx, self.state)
        self.enterRule(localctx, 0, self.RULE_start)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 28
            self.match(IUPACParser.HASH)
            self.state = 35
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==26:
                self.state = 29
                self.match(IUPACParser.LBRACK)
                self.state = 30
                self.branch()
                self.state = 31
                self.match(IUPACParser.RBRACK)
                self.state = 37
                self._errHandler.sync(self)
                _la = self._input.LA(1)

            self.state = 38
            self.begin()
            self.state = 39
            self.match(IUPACParser.HASH)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class BeginContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def branch(self):
            return self.getTypedRuleContext(IUPACParser.BranchContext,0)


        def deriv(self):
            return self.getTypedRuleContext(IUPACParser.DerivContext,0)


        def SPACE(self):
            return self.getToken(IUPACParser.SPACE, 0)

        def TYPE(self):
            return self.getToken(IUPACParser.TYPE, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_begin

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterBegin" ):
                listener.enterBegin(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitBegin" ):
                listener.exitBegin(self)




    def begin(self):

        localctx = IUPACParser.BeginContext(self, self._ctx, self.state)
        self.enterRule(localctx, 2, self.RULE_begin)
        try:
            self.state = 54
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,1,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 41
                self.branch()
                self.state = 42
                self.deriv()
                self.state = 43
                self.match(IUPACParser.SPACE)
                self.state = 44
                self.match(IUPACParser.TYPE)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 46
                self.branch()
                self.state = 47
                self.deriv()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 49
                self.deriv()
                self.state = 50
                self.match(IUPACParser.SPACE)
                self.state = 51
                self.match(IUPACParser.TYPE)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 53
                self.deriv()
                pass


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

        def deriv(self):
            return self.getTypedRuleContext(IUPACParser.DerivContext,0)


        def con(self):
            return self.getTypedRuleContext(IUPACParser.ConContext,0)


        def branch(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(IUPACParser.BranchContext)
            else:
                return self.getTypedRuleContext(IUPACParser.BranchContext,i)


        def LBRACE(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.LBRACE)
            else:
                return self.getToken(IUPACParser.LBRACE, i)

        def RBRACE(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.RBRACE)
            else:
                return self.getToken(IUPACParser.RBRACE, i)

        def getRuleIndex(self):
            return IUPACParser.RULE_branch

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterBranch" ):
                listener.enterBranch(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitBranch" ):
                listener.exitBranch(self)




    def branch(self):

        localctx = IUPACParser.BranchContext(self, self._ctx, self.state)
        self.enterRule(localctx, 4, self.RULE_branch)
        try:
            self.state = 97
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,2,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 56
                self.deriv()
                self.state = 57
                self.con()
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 59
                self.deriv()
                self.state = 60
                self.con()
                self.state = 61
                self.branch()
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 63
                self.match(IUPACParser.LBRACE)
                self.state = 64
                self.branch()
                self.state = 65
                self.match(IUPACParser.RBRACE)
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 67
                self.deriv()
                self.state = 68
                self.con()
                self.state = 69
                self.match(IUPACParser.LBRACE)
                self.state = 70
                self.branch()
                self.state = 71
                self.match(IUPACParser.RBRACE)
                self.state = 72
                self.branch()
                pass

            elif la_ == 5:
                self.enterOuterAlt(localctx, 5)
                self.state = 74
                self.deriv()
                self.state = 75
                self.con()
                self.state = 76
                self.match(IUPACParser.LBRACE)
                self.state = 77
                self.branch()
                self.state = 78
                self.match(IUPACParser.RBRACE)
                self.state = 79
                self.match(IUPACParser.LBRACE)
                self.state = 80
                self.branch()
                self.state = 81
                self.match(IUPACParser.RBRACE)
                self.state = 82
                self.branch()
                pass

            elif la_ == 6:
                self.enterOuterAlt(localctx, 6)
                self.state = 84
                self.deriv()
                self.state = 85
                self.con()
                self.state = 86
                self.match(IUPACParser.LBRACE)
                self.state = 87
                self.branch()
                self.state = 88
                self.match(IUPACParser.RBRACE)
                self.state = 89
                self.match(IUPACParser.LBRACE)
                self.state = 90
                self.branch()
                self.state = 91
                self.match(IUPACParser.RBRACE)
                self.state = 92
                self.match(IUPACParser.LBRACE)
                self.state = 93
                self.branch()
                self.state = 94
                self.match(IUPACParser.RBRACE)
                self.state = 95
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

        def modi(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(IUPACParser.ModiContext)
            else:
                return self.getTypedRuleContext(IUPACParser.ModiContext,i)


        def saci(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(IUPACParser.SaciContext)
            else:
                return self.getTypedRuleContext(IUPACParser.SaciContext,i)


        def RING(self):
            return self.getToken(IUPACParser.RING, 0)

        def TYPE(self):
            return self.getToken(IUPACParser.TYPE, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_deriv

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterDeriv" ):
                listener.enterDeriv(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitDeriv" ):
                listener.exitDeriv(self)




    def deriv(self):

        localctx = IUPACParser.DerivContext(self, self._ctx, self.state)
        self.enterRule(localctx, 6, self.RULE_deriv)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 102
            self._errHandler.sync(self)
            _alt = self._interp.adaptivePredict(self._input,3,self._ctx)
            while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                if _alt==1:
                    self.state = 99
                    self.modi() 
                self.state = 104
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,3,self._ctx)

            self.state = 106 
            self._errHandler.sync(self)
            _alt = 1
            while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                if _alt == 1:
                    self.state = 105
                    self.saci()

                else:
                    raise NoViableAltException(self)
                self.state = 108 
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,4,self._ctx)

            self.state = 113
            self._errHandler.sync(self)
            _alt = self._interp.adaptivePredict(self._input,5,self._ctx)
            while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                if _alt==1:
                    self.state = 110
                    self.modi() 
                self.state = 115
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,5,self._ctx)

            self.state = 117
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            if _la==15:
                self.state = 116
                self.match(IUPACParser.RING)


            self.state = 122
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while (((_la) & ~0x3f) == 0 and ((1 << _la) & 2174448) != 0):
                self.state = 119
                self.modi()
                self.state = 124
                self._errHandler.sync(self)
                _la = self._input.LA(1)

            self.state = 126
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,8,self._ctx)
            if la_ == 1:
                self.state = 125
                self.match(IUPACParser.TYPE)


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
            return self.getToken(IUPACParser.COUNT, 0)

        def SAC(self):
            return self.getToken(IUPACParser.SAC, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_saci

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSaci" ):
                listener.enterSaci(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSaci" ):
                listener.exitSaci(self)




    def saci(self):

        localctx = IUPACParser.SaciContext(self, self._ctx, self.state)
        self.enterRule(localctx, 8, self.RULE_saci)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 128
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

        def LPAR(self):
            return self.getToken(IUPACParser.LPAR, 0)

        def typi(self):
            return self.getTypedRuleContext(IUPACParser.TypiContext,0)


        def NUM(self):
            return self.getToken(IUPACParser.NUM, 0)

        def DASH(self):
            return self.getToken(IUPACParser.DASH, 0)

        def qnum(self):
            return self.getTypedRuleContext(IUPACParser.QnumContext,0)


        def RPAR(self):
            return self.getToken(IUPACParser.RPAR, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_con

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterCon" ):
                listener.enterCon(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitCon" ):
                listener.exitCon(self)




    def con(self):

        localctx = IUPACParser.ConContext(self, self._ctx, self.state)
        self.enterRule(localctx, 10, self.RULE_con)
        try:
            self.state = 151
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,9,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 130
                self.match(IUPACParser.LPAR)
                self.state = 131
                self.typi()
                self.state = 132
                self.match(IUPACParser.NUM)
                self.state = 133
                self.match(IUPACParser.DASH)
                self.state = 134
                self.qnum()
                self.state = 135
                self.match(IUPACParser.RPAR)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 137
                self.match(IUPACParser.LPAR)
                self.state = 138
                self.match(IUPACParser.NUM)
                self.state = 139
                self.match(IUPACParser.DASH)
                self.state = 140
                self.qnum()
                self.state = 141
                self.match(IUPACParser.RPAR)
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 143
                self.typi()
                self.state = 144
                self.match(IUPACParser.NUM)
                self.state = 145
                self.match(IUPACParser.DASH)
                self.state = 146
                self.qnum()
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 148
                self.typi()
                self.state = 149
                self.match(IUPACParser.NUM)
                pass


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
            return self.getToken(IUPACParser.EQ, 0)

        def LBRACK(self):
            return self.getToken(IUPACParser.LBRACK, 0)

        def NUM(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.NUM)
            else:
                return self.getToken(IUPACParser.NUM, i)

        def RBRACK(self):
            return self.getToken(IUPACParser.RBRACK, 0)

        def ct(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(IUPACParser.CtContext)
            else:
                return self.getTypedRuleContext(IUPACParser.CtContext,i)


        def COMMA(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.COMMA)
            else:
                return self.getToken(IUPACParser.COMMA, i)

        def C(self):
            return self.getToken(IUPACParser.C, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_add

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterAdd" ):
                listener.enterAdd(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitAdd" ):
                listener.exitAdd(self)




    def add(self):

        localctx = IUPACParser.AddContext(self, self._ctx, self.state)
        self.enterRule(localctx, 12, self.RULE_add)
        self._la = 0 # Token type
        try:
            self.state = 181
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [34]:
                self.enterOuterAlt(localctx, 1)
                self.state = 153
                self.match(IUPACParser.EQ)
                self.state = 154
                self.match(IUPACParser.LBRACK)
                self.state = 156
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==30 or _la==33:
                    self.state = 155
                    self.ct()


                self.state = 158
                self.match(IUPACParser.NUM)
                self.state = 166
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==18:
                    self.state = 159
                    self.match(IUPACParser.COMMA)
                    self.state = 161
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)
                    if _la==30 or _la==33:
                        self.state = 160
                        self.ct()


                    self.state = 163
                    self.match(IUPACParser.NUM)
                    self.state = 168
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                self.state = 169
                self.match(IUPACParser.RBRACK)
                pass
            elif token in [30]:
                self.enterOuterAlt(localctx, 2)
                self.state = 170
                self.match(IUPACParser.C)
                self.state = 171
                self.match(IUPACParser.LBRACK)
                self.state = 172
                self.match(IUPACParser.NUM)
                self.state = 177
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==18:
                    self.state = 173
                    self.match(IUPACParser.COMMA)
                    self.state = 174
                    self.match(IUPACParser.NUM)
                    self.state = 179
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                self.state = 180
                self.match(IUPACParser.RBRACK)
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
            return self.getToken(IUPACParser.COUNT, 0)

        def bridge(self):
            return self.getTypedRuleContext(IUPACParser.BridgeContext,0)


        def FG(self):
            return self.getToken(IUPACParser.FG, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_fgi

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterFgi" ):
                listener.enterFgi(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitFgi" ):
                listener.exitFgi(self)




    def fgi(self):

        localctx = IUPACParser.FgiContext(self, self._ctx, self.state)
        self.enterRule(localctx, 14, self.RULE_fgi)
        try:
            self.state = 186
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [13]:
                self.enterOuterAlt(localctx, 1)
                self.state = 183
                self.match(IUPACParser.COUNT)
                pass
            elif token in [4, 5, 6, 7]:
                self.enterOuterAlt(localctx, 2)
                self.state = 184
                self.bridge()
                pass
            elif token in [8]:
                self.enterOuterAlt(localctx, 3)
                self.state = 185
                self.match(IUPACParser.FG)
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
            return self.getToken(IUPACParser.CARBON, 0)

        def NUM(self):
            return self.getToken(IUPACParser.NUM, 0)

        def add(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(IUPACParser.AddContext)
            else:
                return self.getTypedRuleContext(IUPACParser.AddContext,i)


        def AI(self):
            return self.getToken(IUPACParser.AI, 0)

        def A(self):
            return self.getToken(IUPACParser.A, 0)

        def I(self):
            return self.getToken(IUPACParser.I, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_carb

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterCarb" ):
                listener.enterCarb(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitCarb" ):
                listener.exitCarb(self)




    def carb(self):

        localctx = IUPACParser.CarbContext(self, self._ctx, self.state)
        self.enterRule(localctx, 16, self.RULE_carb)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 189
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            if (((_la) & ~0x3f) == 0 and ((1 << _la) & 35165044736) != 0):
                self.state = 188
                _la = self._input.LA(1)
                if not((((_la) & ~0x3f) == 0 and ((1 << _la) & 35165044736) != 0)):
                    self._errHandler.recoverInline(self)
                else:
                    self._errHandler.reportMatch(self)
                    self.consume()


            self.state = 191
            self.match(IUPACParser.CARBON)
            self.state = 192
            self.match(IUPACParser.NUM)
            self.state = 196
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==30 or _la==34:
                self.state = 193
                self.add()
                self.state = 198
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
                return self.getTokens(IUPACParser.NUM)
            else:
                return self.getToken(IUPACParser.NUM, i)

        def DASH(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.DASH)
            else:
                return self.getToken(IUPACParser.DASH, i)

        def ANHYDRO(self):
            return self.getToken(IUPACParser.ANHYDRO, 0)

        def COMMA(self):
            return self.getToken(IUPACParser.COMMA, 0)

        def bridge(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(IUPACParser.BridgeContext)
            else:
                return self.getTypedRuleContext(IUPACParser.BridgeContext,i)


        def fgi(self):
            return self.getTypedRuleContext(IUPACParser.FgiContext,0)


        def D(self):
            return self.getToken(IUPACParser.D, 0)

        def E(self):
            return self.getToken(IUPACParser.E, 0)

        def U(self):
            return self.getToken(IUPACParser.U, 0)

        def carb(self):
            return self.getTypedRuleContext(IUPACParser.CarbContext,0)


        def HEAD(self):
            return self.getToken(IUPACParser.HEAD, 0)

        def HEADD(self):
            return self.getToken(IUPACParser.HEADD, 0)

        def END(self):
            return self.getToken(IUPACParser.END, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_modi

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterModi" ):
                listener.enterModi(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitModi" ):
                listener.exitModi(self)




    def modi(self):

        localctx = IUPACParser.ModiContext(self, self._ctx, self.state)
        self.enterRule(localctx, 18, self.RULE_modi)
        self._la = 0 # Token type
        try:
            self.state = 243
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,24,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 201
                self._errHandler.sync(self)
                la_ = self._interp.adaptivePredict(self._input,18,self._ctx)
                if la_ == 1:
                    self.state = 199
                    self.match(IUPACParser.NUM)
                    self.state = 200
                    self.match(IUPACParser.COMMA)


                self.state = 203
                self.match(IUPACParser.NUM)
                self.state = 204
                self.match(IUPACParser.DASH)
                self.state = 205
                self.match(IUPACParser.ANHYDRO)
                self.state = 206
                self.match(IUPACParser.DASH)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 207
                self.match(IUPACParser.NUM)
                self.state = 208
                self.match(IUPACParser.DASH)
                self.state = 209
                self.bridge()
                self.state = 210
                self.match(IUPACParser.DASH)
                self.state = 211
                self.fgi()
                self.state = 212
                self.match(IUPACParser.DASH)
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 215
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==21:
                    self.state = 214
                    self.match(IUPACParser.DASH)


                self.state = 218
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==16:
                    self.state = 217
                    self.match(IUPACParser.NUM)


                self.state = 223
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,21,self._ctx)
                while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                    if _alt==1:
                        self.state = 220
                        self.bridge() 
                    self.state = 225
                    self._errHandler.sync(self)
                    _alt = self._interp.adaptivePredict(self._input,21,self._ctx)

                self.state = 226
                self.fgi()
                pass

            elif la_ == 4:
                self.enterOuterAlt(localctx, 4)
                self.state = 227
                self.match(IUPACParser.NUM)
                self.state = 236
                self._errHandler.sync(self)
                token = self._input.LA(1)
                if token in [18, 31]:
                    self.state = 230
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)
                    if _la==18:
                        self.state = 228
                        self.match(IUPACParser.COMMA)
                        self.state = 229
                        self.match(IUPACParser.NUM)


                    self.state = 232
                    self.match(IUPACParser.D)
                    pass
                elif token in [32]:
                    self.state = 233
                    self.match(IUPACParser.E)
                    pass
                elif token in [36]:
                    self.state = 234
                    self.match(IUPACParser.U)
                    pass
                elif token in [4, 28, 29, 35]:
                    self.state = 235
                    self.carb()
                    pass
                else:
                    raise NoViableAltException(self)

                pass

            elif la_ == 5:
                self.enterOuterAlt(localctx, 5)
                self.state = 238
                self.match(IUPACParser.HEAD)
                pass

            elif la_ == 6:
                self.enterOuterAlt(localctx, 6)
                self.state = 239
                self.match(IUPACParser.HEADD)
                self.state = 240
                self.match(IUPACParser.DASH)
                pass

            elif la_ == 7:
                self.enterOuterAlt(localctx, 7)
                self.state = 241
                self.match(IUPACParser.DASH)
                self.state = 242
                self.match(IUPACParser.END)
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
            return self.getToken(IUPACParser.QMARK, 0)

        def NUM(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.NUM)
            else:
                return self.getToken(IUPACParser.NUM, i)

        def SLASH(self, i:int=None):
            if i is None:
                return self.getTokens(IUPACParser.SLASH)
            else:
                return self.getToken(IUPACParser.SLASH, i)

        def getRuleIndex(self):
            return IUPACParser.RULE_qnum

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterQnum" ):
                listener.enterQnum(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitQnum" ):
                listener.exitQnum(self)




    def qnum(self):

        localctx = IUPACParser.QnumContext(self, self._ctx, self.state)
        self.enterRule(localctx, 20, self.RULE_qnum)
        self._la = 0 # Token type
        try:
            self.state = 254
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [41]:
                self.enterOuterAlt(localctx, 1)
                self.state = 245
                self.match(IUPACParser.QMARK)
                pass
            elif token in [16]:
                self.enterOuterAlt(localctx, 2)
                self.state = 246
                self.match(IUPACParser.NUM)
                self.state = 251
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==37:
                    self.state = 247
                    self.match(IUPACParser.SLASH)
                    self.state = 248
                    self.match(IUPACParser.NUM)
                    self.state = 253
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
            return self.getToken(IUPACParser.TYPE, 0)

        def QMARK(self):
            return self.getToken(IUPACParser.QMARK, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_typi

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterTypi" ):
                listener.enterTypi(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitTypi" ):
                listener.exitTypi(self)




    def typi(self):

        localctx = IUPACParser.TypiContext(self, self._ctx, self.state)
        self.enterRule(localctx, 22, self.RULE_typi)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 256
            _la = self._input.LA(1)
            if not(_la==14 or _la==41):
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
            return self.getToken(IUPACParser.CARBON, 0)

        def NITROGEN(self):
            return self.getToken(IUPACParser.NITROGEN, 0)

        def OXYGEN(self):
            return self.getToken(IUPACParser.OXYGEN, 0)

        def PHOSPHOR(self):
            return self.getToken(IUPACParser.PHOSPHOR, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_bridge

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterBridge" ):
                listener.enterBridge(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitBridge" ):
                listener.exitBridge(self)




    def bridge(self):

        localctx = IUPACParser.BridgeContext(self, self._ctx, self.state)
        self.enterRule(localctx, 24, self.RULE_bridge)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 258
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
            return self.getToken(IUPACParser.C, 0)

        def T(self):
            return self.getToken(IUPACParser.T, 0)

        def getRuleIndex(self):
            return IUPACParser.RULE_ct

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterCt" ):
                listener.enterCt(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitCt" ):
                listener.exitCt(self)




    def ct(self):

        localctx = IUPACParser.CtContext(self, self._ctx, self.state)
        self.enterRule(localctx, 26, self.RULE_ct)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 260
            _la = self._input.LA(1)
            if not(_la==30 or _la==33):
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





