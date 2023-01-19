from setuptools import setup

setup(
    name="GlyLES",
    version="0.5.6",
    description="Convert the IUPAC into a SMILES representation",
    author="Roman Joeres",
    packages=["glyles"],
    scripts=["glyles/glyles"]
)