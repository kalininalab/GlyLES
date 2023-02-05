from antlr4.error.ErrorListener import ErrorListener
from glyles.glycans.utils import ParseError
import sys


class GlyLESErrorListener(ErrorListener):

    def __init__(self):
        super(GlyLESErrorListener, self).__init__()

    def syntaxError(self, recognizer, offendingSymbol, line, column, msg, e):
        error_msg = f"line {line}:{column} {msg}"
        # Printing to stderr kept to mimic the old behaviour.
        # Maybe not a good idea for a library to pollute the stderr when the exceptions are handled?
        print(error_msg, file=sys.stderr)
        raise ParseError(f"Glycan cannot be parsed:\n{error_msg}")
