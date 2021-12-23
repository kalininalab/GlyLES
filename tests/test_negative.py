import sys

from glyles.converter import convert


class TestNegatives:
    def test_invalid_iupac_1(self):
        log = []

        class writer(object):
            def write(self, data):
                log.append(data)

        sys.stderr = writer()

        output = convert("Glca", returning=True)

        assert len(output) == 1
        assert len(output[0]) == 2
        assert output[0][0] == "Glca"
        assert output[0][1] == ""
        assert len(log) == 2
        assert log[0] == "Glycan Glca is not parsable due to unknown monomers."
        assert log[1] == "\n"
