# Generated from Glycan.g4 by ANTLR 4.9
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys



def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\2\b")
        buf.write("/\b\1\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t\7")
        buf.write("\3\2\3\2\3\3\3\3\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3")
        buf.write("\4\3\4\3\4\3\4\3\4\3\4\5\4#\n\4\3\5\3\5\3\5\3\5\3\5\3")
        buf.write("\5\3\5\3\6\3\6\3\7\3\7\2\2\b\3\3\5\4\7\5\t\6\13\7\r\b")
        buf.write("\3\2\2\2\62\2\3\3\2\2\2\2\5\3\2\2\2\2\7\3\2\2\2\2\t\3")
        buf.write("\2\2\2\2\13\3\2\2\2\2\r\3\2\2\2\3\17\3\2\2\2\5\21\3\2")
        buf.write("\2\2\7\"\3\2\2\2\t$\3\2\2\2\13+\3\2\2\2\r-\3\2\2\2\17")
        buf.write("\20\7]\2\2\20\4\3\2\2\2\21\22\7_\2\2\22\6\3\2\2\2\23\24")
        buf.write("\7I\2\2\24\25\7n\2\2\25#\7e\2\2\26\27\7I\2\2\27\30\7n")
        buf.write("\2\2\30#\7w\2\2\31\32\7H\2\2\32\33\7t\2\2\33#\7w\2\2\34")
        buf.write("\35\7O\2\2\35\36\7c\2\2\36#\7p\2\2\37 \7I\2\2 !\7c\2\2")
        buf.write("!#\7n\2\2\"\23\3\2\2\2\"\26\3\2\2\2\"\31\3\2\2\2\"\34")
        buf.write("\3\2\2\2\"\37\3\2\2\2#\b\3\2\2\2$%\7*\2\2%&\5\13\6\2&")
        buf.write("\'\5\r\7\2\'(\7/\2\2()\5\r\7\2)*\7+\2\2*\n\3\2\2\2+,\4")
        buf.write("cd\2,\f\3\2\2\2-.\4\63;\2.\16\3\2\2\2\4\2\"\2")
        return buf.getvalue()


class GlycanLexer(Lexer):

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    T__0 = 1
    T__1 = 2
    SAC = 3
    CON = 4
    TYPE = 5
    NUM = 6

    channelNames = [ u"DEFAULT_TOKEN_CHANNEL", u"HIDDEN" ]

    modeNames = [ "DEFAULT_MODE" ]

    literalNames = [ "<INVALID>",
            "'['", "']'" ]

    symbolicNames = [ "<INVALID>",
            "SAC", "CON", "TYPE", "NUM" ]

    ruleNames = [ "T__0", "T__1", "SAC", "CON", "TYPE", "NUM" ]

    grammarFileName = "Glycan.g4"

    def __init__(self, input=None, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.9")
        self._interp = LexerATNSimulator(self, self.atn, self.decisionsToDFA, PredictionContextCache())
        self._actions = None
        self._predicates = None


