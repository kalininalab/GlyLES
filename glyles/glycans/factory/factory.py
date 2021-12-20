from glyles.glycans.factory.factory_d import DerivativesFactory
from glyles.glycans.factory.factory_f import FuranoseFactory
from glyles.glycans.factory.factory_p import PyranoseFactory


class MonomerFactory:
    def __init__(self):
        self.pyranoses = PyranoseFactory()
        self.furanoses = FuranoseFactory()
        self.derivatives = DerivativesFactory()

        self.keys = set(self.pyranoses.keys()).union(set(self.derivatives.keys()))

    def __contains__(self, item):
        return item.upper() in self.keys

    def keys(self):
        return self.keys

    def monomers(self):
        return set(x.split("_")[-1] for x in self.keys)

    def __getitem__(self, item):
        furanose = False
        if item[-1] == "p":
            item = item[:-1]
        if item[-1] == "f":
            item = item[:-1]
            furanose = True

        if furanose and item in self.furanoses:
            return self.furanoses[item]
        if not furanose and item in self.pyranoses:
            return self.pyranoses[item]
        return self.derivatives[item]
