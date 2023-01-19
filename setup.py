from setuptools import setup

setup(
    name="glyles",
    version="0.5.5",
    description="Convert the IUPAC into a SMILES representation",
    author="Roman Joeres",
    packages=["GlyLES", "converter"],
    scripts=["glyles/glyles"]
)