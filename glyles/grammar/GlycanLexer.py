# Generated from Glycan.g4 by ANTLR 4.9
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys



def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\2\n")
        buf.write("\62\b\1\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t")
        buf.write("\7\4\b\t\b\4\t\t\t\3\2\3\2\3\3\3\3\3\4\3\4\3\5\3\5\3\6")
        buf.write("\3\6\3\7\3\7\3\7\3\7\3\7\3\7\3\7\3\7\3\7\3\7\3\7\3\7\3")
        buf.write("\7\3\7\3\7\5\7-\n\7\3\b\3\b\3\t\3\t\2\2\n\3\3\5\4\7\5")
        buf.write("\t\6\13\7\r\b\17\t\21\n\3\2\2\2\65\2\3\3\2\2\2\2\5\3\2")
        buf.write("\2\2\2\7\3\2\2\2\2\t\3\2\2\2\2\13\3\2\2\2\2\r\3\2\2\2")
        buf.write("\2\17\3\2\2\2\2\21\3\2\2\2\3\23\3\2\2\2\5\25\3\2\2\2\7")
        buf.write("\27\3\2\2\2\t\31\3\2\2\2\13\33\3\2\2\2\r,\3\2\2\2\17.")
        buf.write("\3\2\2\2\21\60\3\2\2\2\23\24\7]\2\2\24\4\3\2\2\2\25\26")
        buf.write("\7_\2\2\26\6\3\2\2\2\27\30\7*\2\2\30\b\3\2\2\2\31\32\7")
        buf.write("/\2\2\32\n\3\2\2\2\33\34\7+\2\2\34\f\3\2\2\2\35\36\7I")
        buf.write("\2\2\36\37\7n\2\2\37-\7e\2\2 !\7I\2\2!\"\7n\2\2\"-\7w")
        buf.write("\2\2#$\7V\2\2$%\7c\2\2%-\7n\2\2&\'\7O\2\2\'(\7c\2\2(-")
        buf.write("\7p\2\2)*\7I\2\2*+\7c\2\2+-\7n\2\2,\35\3\2\2\2, \3\2\2")
        buf.write("\2,#\3\2\2\2,&\3\2\2\2,)\3\2\2\2-\16\3\2\2\2./\4cd\2/")
        buf.write("\20\3\2\2\2\60\61\4\63;\2\61\22\3\2\2\2\4\2,\2")
        return buf.getvalue()


class GlycanLexer(Lexer):

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    T__0 = 1
    T__1 = 2
    T__2 = 3
    T__3 = 4
    T__4 = 5
    SAC = 6
    TYPE = 7
    NUM = 8

    channelNames = [ u"DEFAULT_TOKEN_CHANNEL", u"HIDDEN" ]

    modeNames = [ "DEFAULT_MODE" ]

    literalNames = [ "<INVALID>",
            "'['", "']'", "'('", "'-'", "')'" ]

    symbolicNames = [ "<INVALID>",
            "SAC", "TYPE", "NUM" ]

    ruleNames = [ "T__0", "T__1", "T__2", "T__3", "T__4", "SAC", "TYPE", 
                  "NUM" ]

    grammarFileName = "Glycan.g4"

    def __init__(self, input=None, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.9")
        self._interp = LexerATNSimulator(self, self.atn, self.decisionsToDFA, PredictionContextCache())
        self._actions = None
        self._predicates = None


