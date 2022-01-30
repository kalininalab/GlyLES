from glyles.glycans.factory.factory_f import FuranoseFactory
from glyles.glycans.factory.factory_p import PyranoseFactory
from glyles.glycans.nx_monomer import NXMonomer
from glyles.glycans.rdkit_monomer import RDKitMonomer
from glyles.glycans.utils import Mode, UnreachableError
from glyles.grammar.GlycanLexer import GlycanLexer


class MonomerFactory:
    """
    Class holding and managing the access and generation of all monomers implemented in this package.
    """

    def __init__(self):
        """
        Initialize this factory by creating instances of all "sub-"factories
        """
        self.pyranoses_fac = PyranoseFactory()
        self.furanoses_fac = FuranoseFactory()

        self.keys = set(self.pyranoses_fac.keys())

    def __contains__(self, item):
        """
        Check if an item is part of this factory and can be returned.

        Args:
            item (str): Code of a monomer to be checked

        Returns:
            True if the item is included in the current version of this package
        """
        return item.upper() in self.keys

    def __getitem__(self, item):
        """
        Get an instance of a monomer from this factory.

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

        if furanose and item in self.furanoses_fac:
            return self.furanoses_fac[item]
        if not furanose and item in self.pyranoses_fac:
            return self.pyranoses_fac[item]
        raise NotImplementedError("Query-monomer is neither in pyranoses nor in furanoses")

    def keys(self):
        """
        Get all monomers that are included in this package by their extended name.

        Returns:
            Set of names for all monomers, their available derivatives and configurations (alpha/beta/undefined)
        """
        return self.keys

    def monomers(self):
        """
        Get the names of all monomers in this package ignoring the alpha/beta conformations.

        Returns:
            List, sorted from long to short, of all monomer names in upper case
        """
        return sorted(set([self[x.split("_")[-1]]["name"] for x in self.keys]), key=lambda x: -len(x))

    def furanoses(self):
        """
        Get the names of all monomers available as furanoses in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all furanose names in upper case
        """
        return sorted(set([self[x.split("_")[-1]]["name"] for x in self.furanoses_fac.keys()]))

    def pyranoses(self):
        """
        Get the names of all monomers available as pyranoses in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all pyranose names in upper case
        """
        return sorted(set([self[x.split("_")[-1]]["name"] for x in self.pyranoses_fac.keys()]))

    def monomer_names(self):
        """
        Get the names of all monomers in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in IUPAC notation
        """
        output = set()
        for item in self.monomers():
            if item in self.pyranoses_fac:
                output.add(self.pyranoses_fac[item]["name"])
        return list(output)

    def pyranose_names(self):
        """
        Get the names of all pyranoses in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in IUPAC notation
        """
        output = set()
        for item in self.pyranoses():
            if item in self.pyranoses_fac:
                output.add(self.pyranoses_fac[item]["name"])
        return list(output)

    def furanose_names(self):
        """
        Get the names of all furanoses in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in IUPAC notation
        """
        output = set()
        for item in self.monomers():
            if item in self.furanoses_fac:
                output.add(self.furanoses_fac[item]["name"])
        return list(output)

    def create(self, recipe, mode=Mode.RDKIT_MODE, config=None):
        """
        Create a monomer from its describing IUPAC string with all added side chains.

        Args:
            recipe (List[Tuple[str, int]]): List of modifications, conformations and the root monomer
            mode (Mode): implementation used to represent monomers
            config (str): configuration if monomer is alpha or beta monomer

        Returns:
            Monomer-Instance containing all modifications given in the input
        """
        # determine the class to use to represent the monomers
        if mode == Mode.NETWORKX_MODE:
            monomer_class = NXMonomer
        elif mode == Mode.RDKIT_MODE:
            monomer_class = RDKitMonomer
        else:
            raise ValueError("Unknown representation mode for monomers!")

        # extract key information from the input, i.e. the type, the configuration and pyranose/furanose
        tmp = list(zip(*recipe))
        name = recipe[tmp[1].index(GlycanLexer.SAC)][0]
        config_index = tmp[1].index(GlycanLexer.TYPE) if GlycanLexer.TYPE in tmp[1] else None
        ring_index = tmp[1].index(GlycanLexer.RING) if GlycanLexer.RING in tmp[1] else None

        # generate the full name that is looked up in the factory
        if config is not None and len(config) > 0:
            name = config + "_" + name
        elif config_index is not None:
            name = recipe[config_index][0] + "_" + name

        # get the monomer from the factory
        try:
            if ring_index is not None and recipe[ring_index][0] == "f":
                monomer = monomer_class(**self.furanoses_fac[name], recipe=recipe)
            else:
                monomer = monomer_class(**self.pyranoses_fac[name], recipe=recipe)
        except KeyError:
            raise UnreachableError("Invalid monomer found, should never get so far.")

        # create the final molecule using the molecule's react-method augmented with the recipe of the molecule
        monomer = monomer.react(*tmp)

        return monomer
