# Generated from Glycan.g4 by ANTLR 4.9
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys



def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\2\t")
        buf.write("q\b\1\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t\7")
        buf.write("\4\b\t\b\3\2\3\2\3\3\3\3\3\4\3\4\3\5\3\5\3\5\3\5\3\5\3")
        buf.write("\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5")
        buf.write("\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3")
        buf.write("\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5")
        buf.write("\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3")
        buf.write("\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\5")
        buf.write("\3\5\3\5\3\5\3\5\5\5e\n\5\3\6\3\6\3\6\3\6\3\6\3\6\3\6")
        buf.write("\3\7\3\7\3\b\3\b\2\2\t\3\3\5\4\7\5\t\6\13\7\r\b\17\t\3")
        buf.write("\2\2\2~\2\3\3\2\2\2\2\5\3\2\2\2\2\7\3\2\2\2\2\t\3\2\2")
        buf.write("\2\2\13\3\2\2\2\2\r\3\2\2\2\2\17\3\2\2\2\3\21\3\2\2\2")
        buf.write("\5\23\3\2\2\2\7\25\3\2\2\2\td\3\2\2\2\13f\3\2\2\2\rm\3")
        buf.write("\2\2\2\17o\3\2\2\2\21\22\7\"\2\2\22\4\3\2\2\2\23\24\7")
        buf.write("]\2\2\24\6\3\2\2\2\25\26\7_\2\2\26\b\3\2\2\2\27\30\7H")
        buf.write("\2\2\30\31\7w\2\2\31e\7e\2\2\32\33\7I\2\2\33\34\7c\2\2")
        buf.write("\34e\7n\2\2\35\36\7I\2\2\36\37\7c\2\2\37 \7n\2\2 !\7\65")
        buf.write("\2\2!e\7U\2\2\"#\7I\2\2#$\7c\2\2$%\7n\2\2%&\78\2\2&e\7")
        buf.write("U\2\2\'(\7I\2\2()\7c\2\2)*\7n\2\2*+\7P\2\2+,\7C\2\2,e")
        buf.write("\7e\2\2-.\7I\2\2./\7c\2\2/\60\7n\2\2\60\61\7P\2\2\61\62")
        buf.write("\7C\2\2\62\63\7e\2\2\63\64\7\66\2\2\64e\7U\2\2\65\66\7")
        buf.write("I\2\2\66\67\7c\2\2\678\7n\2\289\7P\2\29:\7C\2\2:;\7e\2")
        buf.write("\2;<\78\2\2<e\7U\2\2=>\7I\2\2>?\7n\2\2?e\7e\2\2@A\7I\2")
        buf.write("\2AB\7n\2\2BC\7e\2\2Ce\7C\2\2DE\7I\2\2EF\7n\2\2FG\7e\2")
        buf.write("\2GH\7P\2\2HI\7C\2\2Ie\7e\2\2JK\7I\2\2KL\7n\2\2LM\7e\2")
        buf.write("\2MN\7P\2\2NO\7C\2\2OP\7e\2\2PQ\78\2\2Qe\7U\2\2RS\7O\2")
        buf.write("\2ST\7c\2\2Te\7p\2\2UV\7P\2\2VW\7g\2\2WX\7w\2\2XY\7\67")
        buf.write("\2\2YZ\7C\2\2Ze\7e\2\2[\\\7P\2\2\\]\7g\2\2]^\7w\2\2^_")
        buf.write("\7\67\2\2_`\7I\2\2`e\7e\2\2ab\7V\2\2bc\7c\2\2ce\7n\2\2")
        buf.write("d\27\3\2\2\2d\32\3\2\2\2d\35\3\2\2\2d\"\3\2\2\2d\'\3\2")
        buf.write("\2\2d-\3\2\2\2d\65\3\2\2\2d=\3\2\2\2d@\3\2\2\2dD\3\2\2")
        buf.write("\2dJ\3\2\2\2dR\3\2\2\2dU\3\2\2\2d[\3\2\2\2da\3\2\2\2e")
        buf.write("\n\3\2\2\2fg\7*\2\2gh\5\r\7\2hi\5\17\b\2ij\7/\2\2jk\5")
        buf.write("\17\b\2kl\7+\2\2l\f\3\2\2\2mn\4cd\2n\16\3\2\2\2op\4\63")
        buf.write(";\2p\20\3\2\2\2\4\2d\2")
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


