import sys
import pytest

from glyles.converter import convert


class TestNegatives:
    @pytest.mark.parametrize("case", ["Glc(c1-2)Glc", "Glc(ab1-2)Glc", "Glc(a11-2)Glc", "Glc(a 1-2)Glc",
                                      "Glc(a1-21)Glc", "Glc{a1-2}Glc", "Glc[a1-2]Glc", "Glc[a1-2)Glc", "Glc(a1-2]Glc",
                                      "Glc(1-2)Glc", "Glc(a-2)Glc", "Glc(a1-)Glc", "Glc(a12)Glc",
                                      "[Man(b1-4)Glc(a1-2)Glc", "Man](b1-4)Glc(a1-2)Glc", "Man[(b1-4)]Glc(a1-2)Glc",
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
