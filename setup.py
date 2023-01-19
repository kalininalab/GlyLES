from setuptools import setup
from glyles import *

setup(
    name="GlyLES",
    version="0.5.6",
    description="Convert the IUPAC into a SMILES representation",
    author="Roman Joeres",
    packages=["GlyLES"],
    scripts=["glyles/glyles"]
)