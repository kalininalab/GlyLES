Welcome to the documentation of GlyLES
======================================

A tool to convert IUPAC representation of Glycans into SMILES representation.
The code is available onf github: https://github.com/kalininalab/glyles
and can be installed with `pip`: https://pypi.org/project/glyles/

GlyLES documentation
====================

Welcome to the documentation of GlyLES!

.. toctree::
    :glob:
    :maxdepth: 1
    :caption: Notes

    notes/installation
    notes/notation

.. toctree::
    :glob:
    :maxdepth: 1
    :caption: Usage

    howto/conversion
    howto/visualization
    howto/filtering
    howto/fgs

.. toctree::
    :glob:
    :maxdepth: 1
    :caption: Source

    source/convert
    source/glycan

Citation

If you used GlyLES in your work or it was helpful for your project, please cite it.

.. code-block::

    @article{joeres2023glyles,
        title={GlyLES: Grammar-based Parsing of Glycans from IUPAC-condensed to SMILES},
        author={Joeres, Roman and Bojar, Daniel and Kalinina, Olga V},
        journal={Journal of Cheminformatics},
        volume={15},
        number={1},
        pages={1--11},
        year={2023},
        publisher={BioMed Central}
    }
