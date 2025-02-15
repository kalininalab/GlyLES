# Change Log

## v1.2.0 - 2025-02-15

### Added

- GWB support
- IUPAC Export for GWB sequences

## v1.1.0 - 2024-12-18

## v1.0.0 - 2024-08-16

## v0.5.11 - 2023-06-09

Relax requirement restrictions to allow for RDKit versions newer then 2021-09-02

## v0.5.10

## v0.5.2-0.5.9 - 2023-01-26/27

Implementation of a CLI for GlyLES including debugging

### Added

- CLI for GlyLES

## v0.5.1 - 2023-01-12

Release for publication in Journal of Cheminformatics 

### Added

- Incorporate for networkx 3.0 use in glycowork

---

## v0.5.0 - 2022-11-09

Release for publication in Journal of Cheminformatics 

### Added

- ?,?-Anhydro configurations of glycan implemented for structural representation
- Open-form representations for "-aric", "-onic", "-ulosonic", and "-ulosaric" glycans
- Example notebooks and documentation on [ReadTheDocs](https://glyles.readthedocs.io/en/latest/).
- 99% of the structures in the glycowork database can be converted (if they are convertible).

---

## v0.4.2 - 2022-09-06

### Fixed

- Bug fixed

---

## v0.4.1 - 2022-09-02

### Fixed

- Patch for full representation mode of glycans
- Logging levels based on python logging module introduced instead of own logging mix

---

## v0.4.0 - 2022-08-25

### Added

- Ternary and Quarternary branching is now possible. Monomers may have up to four children in the tree-like structure 
of a glycan.
- More tests for the new functionality
- Support for two new IUPAC representations (full and simplified according to 
[SNFG](https://www.ncbi.nlm.nih.gov/glycans/snfg.html)) beside the current IUPAC-condensed
- Now supporting most of the functional groups in 
  [CSDB](http://csdb.glycoscience.ru/snfgedit/snfgedit.html?expert=1&destination=structure) and 
  [glycowork](https://pypi.org/project/glycowork/). 

### Changed

- Monomers only in RDKit representation possible. NetworkX-implementation has been removed
- Adding of functional groups completely reworked
- Comments and code style improved, also a bit of a cleanup
- Default output-type is ``return`` instead of commandline prints

---

## v0.3.1 - 2022-05-18

### Fixed

- Bug fixed in ANTLR4 package version requirement

---

## v0.3.0 - 2022-02-10

Release for 17th German Conference on Cheminformatics and EuroSAMPL Satellite Workshop
([Link](https://veranstaltungen.gdch.de/tms/frontend/index.cfm?l=10916&sp_id=2))

### Added

- Length restrictions for glycans lifted, now also glycan trees of depth > 9 parsable

### Changed

- Improved descriptions in READMEs

---

## v0.2.1 - 2022-02-10

### Added

- Detailed description
- Two more tests for the length of branches
- Changelog

### Changed

- Improved descriptions in README

---

## v0.2.0 - 2022-01-26

### Added

- More tests
- More parsable and translatable modifications

### Changed

- Representation of Neuraminic Acid

### Fixed

- Bug fixes regarding C1-atom detection

---

## v0.1.2 - 2022-01-16

### Added

- More tests

### Fixed

- Bugs with modifications fixed


---

## v0.1.1 - 2022-01-10

### Added

- More tests
- Bugs with modifications reported in README

---

## v0.1.0 - 2022-01-03

### Added

- Distinction in pyranoses and furanoses
- More tests

---

## v0.0.4 - 2021-12-15

### Fixed

- Description in README

---

## v0.0.3 - 2021-12-13

### Fixed

- README made parsable for pypi

---

## v0.0.2 - 2021-12-13

### Changed

- Code style and README improvements

---

## v0.0.1 - 2021-12-13

Initial release

### Added

- Code base to convert SMILES from IUPAC notation to SMILES notation
