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


@pytest.mark.parametrize("gwb", [
    "redEnd--??1D-Glc,p--4b1D-Gal,p(--3a2D-NeuGc,p@270)--4b1D-GalNAc,p$MONO,Und,-H,0,redEnd"
    "redEnd--??1D-Glc,p--4b1D-Gal,p@360(--4b1D-GalNAc,p)--3a2D-NeuGc,p@360$MONO,perMe,Na,0,redEnd"
])
def test_gwb_positioning(gwb):
    assert GWBGlycan(gwb).to_iupac(slim=True) == "NeuGc(a2-3)[GalNAc(b1-4)]Gal(b1-4)Glc"


@pytest.mark.parametrize("data", iupac_data)
def test_gwb_iupac(data):
    gwb, iupac = data
    assert GWBGlycan(gwb).to_iupac(slim=True) == iupac


@pytest.mark.fuzzy
@pytest.mark.parametrize("gwb", fuzzy_data)
def test_gwb_fuzzy(gwb):
    if any(x in gwb for x in {"[", "]"}):  # skip repeats
        pytest.skip()
    iupac = GWBGlycan(gwb).to_iupac()
    assert iupac is not None
    assert len(iupac) >= 3
