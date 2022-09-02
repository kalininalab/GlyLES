from glyles.glycans.factory.factory_f import FuranoseFactory
from glyles.glycans.factory.factory_o import c1_finder, OpenFactory
from glyles.glycans.factory.factory_p import PyranoseFactory
from glyles.glycans.monomer import Monomer
from glyles.glycans.utils import Config, Enantiomer, Lactole
from glyles.grammar.GlycanLexer import GlycanLexer


class MonomerFactory:
    """
    Class holding and managing the access and generation of all monomers implemented in this package.
    """

    def __init__(self):
        """
        Initialize this factory by creating instances of all "sub-"factories
        """
        self.pyranose_fac = PyranoseFactory()
        self.furanose_fac = FuranoseFactory()
        self.open_fac = OpenFactory()

        self.keys = set(self.pyranose_fac.keys()).union(self.furanose_fac.keys())

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
        furanose = item in ["ERY", "THRE", "RUL", "XUL", "XLU", "ACE"]
        if item[-1] == "p" and not item.endswith("Hep"):
            item = item[:-1]
        if item[-1] == "f":
            item = item[:-1]
            furanose = True

        if furanose and item in self.furanose_fac:
            return self.furanose_fac[item]
        if not furanose and item in self.pyranose_fac:
            return self.pyranose_fac[item]
        if item in self.open_fac:
            return self.open_fac[item]
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
        return sorted(set([self[x.split("_")[-1]]["name"] for x in self.furanose_fac.keys()]))

    def pyranoses(self):
        """
        Get the names of all monomers available as pyranoses in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all pyranose names in upper case
        """
        return sorted(set([self[x.split("_")[-1]]["name"] for x in self.pyranose_fac.keys()]))

    def monomer_names(self):
        """
        Get the names of all monomers in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in IUPAC notation
        """
        output = set()
        for item in self.monomers():
            if item in self.pyranose_fac:
                output.add(self.pyranose_fac[item]["name"])
        return list(output)

    def pyranose_names(self):
        """
        Get the names of all pyranoses in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in IUPAC notation
        """
        output = set()
        for item in self.pyranoses():
            if item in self.pyranose_fac:
                output.add(self.pyranose_fac[item]["name"])
        return list(output)

    def furanose_names(self):
        """
        Get the names of all furanoses in this package ignoring the alpha/beta conformations

        Returns:
            List, sorted from long to short, of all monomer names in IUPAC notation
        """
        output = set()
        for item in self.monomers():
            if item in self.furanose_fac:
                output.add(self.furanose_fac[item]["name"])
        return list(output)

    @staticmethod
    def unknown_monomer(name):
        """
        Generate a dummy monosaccharide containing no information except the name as the structure and all its
        properties are unknown.
        Args:
            name (str): Name of the monosaccharide with structure and properties are unknown

        Returns:
            Dummy information about an unknown monomer
        """
        return {"name": name, "config": Config.UNDEF, "isomer": Enantiomer.U, "lactole": Lactole.UNKNOWN,
                "smiles": ""}

    @staticmethod
    def succinic_acid():
        """
        This method stores the information about succinic acid in a dictionary.

        Returns:
            Information about succinic acid
        """
        return {"name": "Suc", "config": Config.UNDEF, "isomer": Enantiomer.U, "lactole": Lactole.OPEN,
                "smiles": "OC[C@H](O)CCO", "c1_find": lambda x: c1_finder(x, "OC[C@H](O)CCO")}

    def create(self, recipe, config=None, tree_only=False):
        """
        Create a monomer from its describing IUPAC string with all added side chains.

        Args:
            recipe (List[Tuple[str, int]]): List of modifications, conformations and the root monomer
            config (str): configuration if monomer is alpha or beta monomer
            tree_only (bool): Flag indicating to only parse the tree of glycans and not the modifications

        Returns:
            Monomer-Instance containing all modifications given in the input
        """

        # extract key information from the input, i.e. the type, the configuration and pyranose/furanose
        tmp = list(zip(*recipe))
        name = recipe[tmp[1].index(GlycanLexer.SAC)][0]
        if name == "Sug":
            name = "Oct"
        config_index = tmp[1].index(GlycanLexer.TYPE) if GlycanLexer.TYPE in tmp[1] else None
        ring_index = tmp[1].index(GlycanLexer.RING) if GlycanLexer.RING in tmp[1] else None

        # generate the full name that is looked up in the factory
        if config is not None and len(config) > 0:
            name = config + "_" + name
        elif config_index is not None:
            name = recipe[config_index][0] + "_" + name

        # get the monomer from the factory
        if name in self.pyranose_fac and (ring_index is None or recipe[ring_index][0] != "f"):
            monomer = Monomer(**self.pyranose_fac[name], recipe=recipe)
        elif name in self.furanose_fac:
            monomer = Monomer(**self.furanose_fac[name], recipe=recipe)
        elif name in self.open_fac:
            monomer = Monomer(**self.open_fac[name], recipe=recipe)
        elif name[-3:].upper() == "SUC":
            monomer = Monomer(**self.succinic_acid(), recipe=recipe)
        else:
            monomer = Monomer(**self.unknown_monomer(name), recipe=recipe)

        full = False
        if not tree_only:
            # create the final molecule using the molecule's react-method augmented with the recipe of the molecule
            monomer, full = monomer.react(*tmp)

        # set full to false in case the monomer is unknown, i.e. neither pyranose nor furanose nor open-form
        full &= (monomer.get_lactole() != Lactole.UNKNOWN)

        return monomer, full
