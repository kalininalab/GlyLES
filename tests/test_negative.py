import pytest

from glyles.converter import convert


class TestNegatives:
    @pytest.mark.parametrize("case", ["Glc(c1-2)Glc", "Glc(ab1-2)Glc", "Glc(a11-2)Glc", "Glc(a 1-2)Glc",
                                      "Glc(a1-21)Glc", "Glc{a1-2}Glc", "Glc[a1-2]Glc", "Glc[a1-2)Glc", "Glc(a1-2]Glc",
                                      "Glc(1-2)Glc", "Glc(a-2)Glc", "Glc(a1-)Glc", "Glc(a12)Glc",
                                      "[Man(b1-4)Glc(a1-2)Glc", "Man](b1-4)Glc(a1-2)Glc", "Man[(b1-4)]Glc(a1-2)Glc",
                                      "[Man(b1]-4)Glc(a1-2)Glc", "[Man(b1-4)Glc](a1-2)Glc", "Man[(b1-4)Glc](a1-2)Glc",
                                      "Man[(b1-4)Glc(a1-2)]Glc", "Man(b1-4)[Glc(a1-2)Glc]", "Man(b1-4)Glc(a1-2)[Glc]",
                                      "Man(b1-4)[Glc(a1-2)Glc]", "[Man(b1-4)Glc(a1-2)Glc]", "Glca(c1-2)Glc",
                                      "Glcb(ab1-2)Glc", "Glca(a11-2)Glc", "Glcb(a 1-2)Glc", "Glca(a1-21)Glc",
                                      "Glcb{a1-2}Glc", "Glca[a1-2]Glc", "Glcb[a1-2)Glc", "Glca(a1-2]Glc",
                                      "Glcb(1-2)Glc", "Glca(a-2)Glc", "Glcb(a1-)Glc", "Glca(a12)Glc", "Glcpa(c1-2)Glc",
                                      "Glcfb(ab1-2)Glc", "Glcfa(a11-2)Glc", "Glcpb(a 1-2)Glc",  "Glcpa(a1-21)Glc",
                                      "Glcfb{a1-2}Glc", "Glcfa[a1-2]Glc", "Glcpb[a1-2)Glc", "Glcpa(a1-2]Glc",
                                      "Glcfb(1-2)Glc", "Glcfa(a-2)Glc", "Glcpb(a1-)Glc", "Glcfa(a12)Glc",
                                      "Glcp(c1-2)Glc", "Glcf(ab1-2)Glc", "Glcf(a11-2)Glc", "Glcp(a 1-2)Glc",
                                      "Glcp(a1-21)Glc", "Glcf{a1-2}Glc", "Glcf[a1-2]Glc", "Glcp[a1-2)Glc",
                                      "Glcp(a1-2]Glc", "Glcf(1-2)Glc", "Glcf(a-2)Glc", "Glcp(a1-)Glc", "Glcf(a12)Glc",
                                      "[Mana(b1-4)Glcb(a1-2)Glc", "Manf](b1-4)Glca(a1-2)Glc",
                                      "Manb[(b1-4)]Glcp(a1-2)Glc", "[Manfb(b1]-4)Glcpa(a1-2)Glc",
                                      "[Manf(b1-4)Glcp](a1-2)Glc", "Manpa[(b1-4)Glcpb](a1-2)Glcf",
                                      "Mana[(b1-4)Glcp(a1-2)]Glcf", "Manf(b1-4)[Glcpb(a1-2)Glcfa]",
                                      "Manf(b1-4)Glca(a1-2)[Glcb]", "Manp(b1-4)[Glcfb(a1-2)Glc]",
                                      "[Manpa(b1-4)Glcf(a1-2)Glcb]"])
    @pytest.mark.parametrize("config", [" a", " b", ""])
    def test_invalid_iupac_full(self, case, config):
        output = convert(case + config, returning=True)

        assert len(output) == 1
        assert len(output[0]) == 2
        assert output[0][0] == case + config
        assert output[0][1] == ""
