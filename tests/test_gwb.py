import re
import pytest
from glyles.glycans.poly.gwb_glycan import GWBGlycan


iupac_data = []
with open("data/gwb_iupac.tsv", "r") as f:
    for line in f.readlines():
        parts = line.strip().split("\t")
        iupac_data.append((parts[0], parts[1]))


fuzzy_data = []
with open("data/gwb_fuzzy.txt", "r") as f:
    for line in f.readlines():
        fuzzy_data.append(line.strip())


@pytest.mark.parametrize("data", iupac_data)
def test_gwb_iupac(data):
    gwb, iupac = data
    assert GWBGlycan(gwb).to_iupac(slim=True) == iupac


@pytest.mark.fuzzy
@pytest.mark.parametrize("gwb", fuzzy_data)
def test_gwb_fuzzy(gwb):
    if any(x in gwb for x in {"dTal", "End--??1P$", "End--??1S$", "[", "]"}):  # skip specific stuff and repeats
        pytest.skip()
    iupac = GWBGlycan(gwb).to_iupac()
    assert iupac is not None
    assert len(iupac) >= 3
