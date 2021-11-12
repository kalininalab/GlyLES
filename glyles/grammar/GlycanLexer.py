# Generated from Glycan.g4 by ANTLR 4.9
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys



def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\2\b")
        buf.write("m\b\1\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t\7")
        buf.write("\3\2\3\2\3\3\3\3\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3")
        buf.write("\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4")
        buf.write("\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3")
        buf.write("\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4")
        buf.write("\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3")
        buf.write("\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4")
        buf.write("\5\4a\n\4\3\5\3\5\3\5\3\5\3\5\3\5\3\5\3\6\3\6\3\7\3\7")
        buf.write("\2\2\b\3\3\5\4\7\5\t\6\13\7\r\b\3\2\2\2z\2\3\3\2\2\2\2")
        buf.write("\5\3\2\2\2\2\7\3\2\2\2\2\t\3\2\2\2\2\13\3\2\2\2\2\r\3")
        buf.write("\2\2\2\3\17\3\2\2\2\5\21\3\2\2\2\7`\3\2\2\2\tb\3\2\2\2")
        buf.write("\13i\3\2\2\2\rk\3\2\2\2\17\20\7]\2\2\20\4\3\2\2\2\21\22")
        buf.write("\7_\2\2\22\6\3\2\2\2\23\24\7H\2\2\24\25\7w\2\2\25a\7e")
        buf.write("\2\2\26\27\7I\2\2\27\30\7c\2\2\30a\7n\2\2\31\32\7I\2\2")
        buf.write("\32\33\7c\2\2\33\34\7n\2\2\34\35\7\65\2\2\35a\7U\2\2\36")
        buf.write("\37\7I\2\2\37 \7c\2\2 !\7n\2\2!\"\78\2\2\"a\7U\2\2#$\7")
        buf.write("I\2\2$%\7c\2\2%&\7n\2\2&\'\7P\2\2\'(\7C\2\2(a\7e\2\2)")
        buf.write("*\7I\2\2*+\7c\2\2+,\7n\2\2,-\7P\2\2-.\7C\2\2./\7e\2\2")
        buf.write("/\60\7\66\2\2\60a\7U\2\2\61\62\7I\2\2\62\63\7c\2\2\63")
        buf.write("\64\7n\2\2\64\65\7P\2\2\65\66\7C\2\2\66\67\7e\2\2\678")
        buf.write("\78\2\28a\7U\2\29:\7I\2\2:;\7n\2\2;a\7e\2\2<=\7I\2\2=")
        buf.write(">\7n\2\2>?\7e\2\2?a\7C\2\2@A\7I\2\2AB\7n\2\2BC\7e\2\2")
        buf.write("CD\7P\2\2DE\7C\2\2Ea\7e\2\2FG\7I\2\2GH\7n\2\2HI\7e\2\2")
        buf.write("IJ\7P\2\2JK\7C\2\2KL\7e\2\2LM\78\2\2Ma\7U\2\2NO\7O\2\2")
        buf.write("OP\7c\2\2Pa\7p\2\2QR\7P\2\2RS\7g\2\2ST\7w\2\2TU\7\67\2")
        buf.write("\2UV\7C\2\2Va\7e\2\2WX\7P\2\2XY\7g\2\2YZ\7w\2\2Z[\7\67")
        buf.write("\2\2[\\\7I\2\2\\a\7e\2\2]^\7V\2\2^_\7c\2\2_a\7n\2\2`\23")
        buf.write("\3\2\2\2`\26\3\2\2\2`\31\3\2\2\2`\36\3\2\2\2`#\3\2\2\2")
        buf.write("`)\3\2\2\2`\61\3\2\2\2`9\3\2\2\2`<\3\2\2\2`@\3\2\2\2`")
        buf.write("F\3\2\2\2`N\3\2\2\2`Q\3\2\2\2`W\3\2\2\2`]\3\2\2\2a\b\3")
        buf.write("\2\2\2bc\7*\2\2cd\5\13\6\2de\5\r\7\2ef\7/\2\2fg\5\r\7")
        buf.write("\2gh\7+\2\2h\n\3\2\2\2ij\4cd\2j\f\3\2\2\2kl\4\63;\2l\16")
        buf.write("\3\2\2\2\4\2`\2")
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


