from glyles.glycans.factory.factory_d import DerivativesFactory
from glyles.glycans.factory.factory_f import FuranoseFactory
from glyles.glycans.factory.factory_p import PyranoseFactory


class MonomerFactory:
    """
    Class holding and managing the access and generation of all monomers implemented in this package
    """

    def __init__(self):
        """
        Initialize this factory by creating instances of all "sub-"factories
        """
        self.pyranoses = PyranoseFactory()
        self.furanoses = FuranoseFactory()
        self.derivatives = DerivativesFactory()

        self.keys = set(self.pyranoses.keys()).union(set(self.derivatives.keys()))

    def __contains__(self, item):
        """
        Check if an item is part of this factory and can be returned

        Args:
            item (str): Code of a monomer to be checked

        Returns:
            True if the item is included in the current version of this package
        """
        return item.upper() in self.keys

    def keys(self):
        """
        Get all monomers that are included in this package by their extended name

        Returns:
            Set of names for all monomers, their available derivatives and configurations (alpha/beta/undefined)
        """
        return self.keys

    def monomers(self):
        """
        Get the names of all monomers in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in upper case
        """
        return sorted([self[x.split("_")[-1]]["name"] for x in self.keys], key=lambda x: -len(x))

    def monomer_names(self):
        """
        Get the names of all monomers in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in IUPAC notation
        """
        output = set()
        for item in self.monomers():
            if item in self.pyranoses:
                output.add(self.pyranoses[item]["name"])
            elif item in self.derivatives:
                output.add(self.derivatives[item]["name"])
        return output

    def __getitem__(self, item):
        """
        Get an instance of a monomer from this factory

        Args:
            item (str): name of the query monomer

        Returns:
            Directory containing all necessary information to initialize a monomer implementation
        """
        furanose = False
        if item[-1] == "p" and not item.endswith("manHep"):
            item = item[:-1]
        if item[-1] == "f":
            item = item[:-1]
            furanose = True

        if furanose and item in self.furanoses:
            return self.furanoses[item]
        if not furanose and item in self.pyranoses:
            return self.pyranoses[item]
        return self.derivatives[item]
