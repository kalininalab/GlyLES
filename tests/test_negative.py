import sys
import pytest

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

    @pytest.mark.parametrize("case", ["Glc(c1-2)Glc", "Glc(ab1-2)Glc", "Glc(a11-2)Glc", "Glc(a 1-2)Glc",
                                      "Glc(a1-21)Glc", "Glc{a1-2}Glc", "Glc[a1-2]Glc", "Glc[a1-2)Glc", "Glc(a1-2]Glc",
                                      "Glc(1-2)Glc", "Glc(a-2)Glc", "Glc(a1-)Glc", "Glc(a12)Glc"])
    @pytest.mark.parametrize("config", [" a", " b", ""])
    def test_invalid_iupac_cons(self, case, config):
        output = convert(case + config, returning=True)

        assert len(output) == 1
        assert len(output[0]) == 2
        assert output[0][0] == case + config
        assert output[0][1] == ""

    @pytest.mark.parametrize("case", ["[Man(b1-4)Glc(a1-2)Glc", "Man](b1-4)Glc(a1-2)Glc", "Man[(b1-4)]Glc(a1-2)Glc",
                                      "[Man(b1]-4)Glc(a1-2)Glc", "[Man(b1-4)Glc](a1-2)Glc", "Man[(b1-4)Glc](a1-2)Glc",
                                      "Man[(b1-4)Glc(a1-2)]Glc", "Man(b1-4)[Glc(a1-2)Glc]", "Man(b1-4)Glc(a1-2)[Glc]",
                                      "Man(b1-4)[Glc(a1-2)Glc]", "[Man(b1-4)Glc(a1-2)Glc]"])
    @pytest.mark.parametrize("config", [" a", " b", ""])
    def test_invalid_iupac_cons(self, case, config):
        output = convert(case + config, returning=True)

        assert len(output) == 1
        assert len(output[0]) == 2
        assert output[0][0] == case + config
        assert output[0][1] == ""
