# Generated from Glycan.g4 by ANTLR 4.9
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys



def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\2\t")
        buf.write("e\b\1\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t\7")
        buf.write("\4\b\t\b\3\2\3\2\3\3\3\3\3\4\3\4\3\5\3\5\3\5\3\5\3\5\3")
        buf.write("\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5")
        buf.write("\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3")
        buf.write("\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5")
        buf.write("\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3")
        buf.write("\5\3\5\3\5\3\5\3\5\3\5\5\5Y\n\5\3\6\3\6\3\6\3\6\3\6\3")
        buf.write("\6\3\6\3\7\3\7\3\b\3\b\2\2\t\3\3\5\4\7\5\t\6\13\7\r\b")
        buf.write("\17\t\3\2\2\2p\2\3\3\2\2\2\2\5\3\2\2\2\2\7\3\2\2\2\2\t")
        buf.write("\3\2\2\2\2\13\3\2\2\2\2\r\3\2\2\2\2\17\3\2\2\2\3\21\3")
        buf.write("\2\2\2\5\23\3\2\2\2\7\25\3\2\2\2\tX\3\2\2\2\13Z\3\2\2")
        buf.write("\2\ra\3\2\2\2\17c\3\2\2\2\21\22\7\"\2\2\22\4\3\2\2\2\23")
        buf.write("\24\7]\2\2\24\6\3\2\2\2\25\26\7_\2\2\26\b\3\2\2\2\27\30")
        buf.write("\7H\2\2\30\31\7w\2\2\31Y\7e\2\2\32\33\7I\2\2\33\34\7c")
        buf.write("\2\2\34Y\7n\2\2\35\36\7I\2\2\36\37\7c\2\2\37 \7n\2\2 ")
        buf.write("!\7\65\2\2!Y\7U\2\2\"#\7I\2\2#$\7c\2\2$%\7n\2\2%&\78\2")
        buf.write("\2&Y\7U\2\2\'(\7I\2\2()\7c\2\2)*\7n\2\2*+\7P\2\2+,\7C")
        buf.write("\2\2,Y\7e\2\2-.\7I\2\2./\7c\2\2/\60\7n\2\2\60\61\7P\2")
        buf.write("\2\61\62\7C\2\2\62\63\7e\2\2\63\64\7\66\2\2\64Y\7U\2\2")
        buf.write("\65\66\7I\2\2\66\67\7c\2\2\678\7n\2\289\7P\2\29:\7C\2")
        buf.write("\2:;\7e\2\2;<\78\2\2<Y\7U\2\2=>\7I\2\2>?\7n\2\2?Y\7e\2")
        buf.write("\2@A\7I\2\2AB\7n\2\2BC\7e\2\2CY\7C\2\2DE\7I\2\2EF\7n\2")
        buf.write("\2FG\7e\2\2GH\7P\2\2HI\7C\2\2IY\7e\2\2JK\7I\2\2KL\7n\2")
        buf.write("\2LM\7e\2\2MN\7P\2\2NO\7C\2\2OP\7e\2\2PQ\78\2\2QY\7U\2")
        buf.write("\2RS\7O\2\2ST\7c\2\2TY\7p\2\2UV\7V\2\2VW\7c\2\2WY\7n\2")
        buf.write("\2X\27\3\2\2\2X\32\3\2\2\2X\35\3\2\2\2X\"\3\2\2\2X\'\3")
        buf.write("\2\2\2X-\3\2\2\2X\65\3\2\2\2X=\3\2\2\2X@\3\2\2\2XD\3\2")
        buf.write("\2\2XJ\3\2\2\2XR\3\2\2\2XU\3\2\2\2Y\n\3\2\2\2Z[\7*\2\2")
        buf.write("[\\\5\r\7\2\\]\5\17\b\2]^\7/\2\2^_\5\17\b\2_`\7+\2\2`")
        buf.write("\f\3\2\2\2ab\4cd\2b\16\3\2\2\2cd\4\63;\2d\20\3\2\2\2\4")
        buf.write("\2X\2")
        return buf.getvalue()


class GlycanLexer(Lexer):

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    T__0 = 1
    T__1 = 2
    T__2 = 3
    SAC = 4
    CON = 5
    TYPE = 6
    NUM = 7

    channelNames = [ u"DEFAULT_TOKEN_CHANNEL", u"HIDDEN" ]

    modeNames = [ "DEFAULT_MODE" ]

    literalNames = [ "<INVALID>",
            "' '", "'['", "']'" ]

    symbolicNames = [ "<INVALID>",
            "SAC", "CON", "TYPE", "NUM" ]

    ruleNames = [ "T__0", "T__1", "T__2", "SAC", "CON", "TYPE", "NUM" ]

    grammarFileName = "Glycan.g4"

    def __init__(self, input=None, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.9")
        self._interp = LexerATNSimulator(self, self.atn, self.decisionsToDFA, PredictionContextCache())
        self._actions = None
        self._predicates = None


