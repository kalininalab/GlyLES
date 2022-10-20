Basic Usage for Glycan Conversion from IUPAC-condensed to SMILES
=====

In the following, we will need those three imports

.. code:: python

    import logging
    from glyles import convert, convert_generator


Input formats
-----

We will with go through different input formats. These include single glycans, lists or tuples, files, and generators of glycans. Combinations of those are possible.


Conversion of a single glycan
^^^^^

The output will be a list with one tuple containing the input IUPAC-condensed string and the SMILES string.

.. code:: python

    convert("Gal")

Output:

.. code:: console

    [('Gal', 'O1C(O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1CO')]

Conversion of a list of glycans
^^^^^

The output will be a list with as many tuples as structures are given. All tuples will hold the input IUPAC-condensed string and the according SMILES string. This works with tuple-input as well.

.. code:: python

    convert(glycan_list=["Gal", "Man", "GalNAc", "NeuAc"])

Output:

.. code:: console

    [('Gal', 'O1C(O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('Man', 'O1C(O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1CO'),
     ('GalNAc', 'O1C(O)[C@H](NC(C)=O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('NeuAc', 'O1[C@@H]([C@H](O)[C@H](O)CO)[C@H](NC(C)=O)[C@@H](O)CC1(O)C(=O)O')]

Conversion of a file containing glycans
^^^^^

Input can also be a file where each file contains one glycan per line.

.. code:: python

    convert(glycan_file="files/general.txt")

Output:

.. code:: console

    [('GalNAc', 'O1C(O)[C@H](NC(C)=O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('GlcNAca', 'O1[C@H](O)[C@H](NC(C)=O)[C@@H](O)[C@H](O)[C@H]1CO'),
     ('GalA', 'O1C(O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1C(=O)O'),
     ('GlcA', 'O1C(O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1C(=O)O'),
     ('GalN', 'O1C(O)[C@H](N)[C@@H](O)[C@@H](O)[C@H]1CO')]

Conversion of a generator of glycans
^^^^^

The output will be the same as if the input is a list of glycans.

.. code:: python

    def gen():
        for glycan in ["Gal", "Man", "GalNAc", "NeuAc"]:
            yield glycan


    convert(glycan_generator=gen())

Output:

.. code:: console

    [('Gal', 'O1C(O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('Man', 'O1C(O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1CO'),
     ('GalNAc', 'O1C(O)[C@H](NC(C)=O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('NeuAc', 'O1[C@@H]([C@H](O)[C@H](O)CO)[C@H](NC(C)=O)[C@@H](O)CC1(O)C(=O)O')]

Output formats
-----

As shown above, we have different input formats, but also the output can be changed. Aside from the list of tuples already shown above, there are other options. Default is to return a list of tuples with the IUPAC-condensed strings and SMILES strings as seen above.

Output as a generator
^^^^^

As output one can select to have a generator which is especially helpful if the input is also a generator and one wants to parallelize jobs. For this specific output, one has to call the convert_generator function. The generator generates tuples of the form as described above.

.. code:: python

    output = convert_generator(glycan_generator=gen())
    print(type(output))
    for g in output:
        print(g)

Output:

.. code:: console

    <class 'generator'>
    ('Gal', 'O1C(O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1CO')
    ('Man', 'O1C(O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1CO')
    ('GalNAc', 'O1C(O)[C@H](NC(C)=O)[C@@H](O)[C@@H](O)[C@H]1CO')
    ('NeuAc', 'O1[C@@H]([C@H](O)[C@H](O)CO)[C@H](NC(C)=O)[C@@H](O)CC1(O)C(=O)O')

Write output to a file
^^^^^

One can write the tuples generated into a file. The file will have the CSV format.

.. code:: python

    convert(glycan_list=["Gal", "Man", "GalNAc", "NeuAc"], output_file="files/example.csv")

Write output to stdout
^^^^^

If file-output is not selected and returning the list is actively rejected, the output will be printed to the standard output. As with the file-output, the output will have CSV format. If some input cannot be converted, the package rises an error and returns a half-empty tuple.

.. code:: python

    convert(glycan_list=["Gal", "Man", "GalNAc", "NeuAc", "Ne{u"], returning=False)

Output:

.. code:: console

    WARNING:root:No output-file specified, results will be printed on stdout.
    ERROR:root:A parsing error occurred with Ne{u: <class 'glyles.glycans.utils.ParseError'>
    Error message: Glycan cannot be parsed:
    line 1:2 token recognition error at: 'e{'
    Gal,O1C(O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1CO
    Man,O1C(O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1CO
    GalNAc,O1C(O)[C@H](NC(C)=O)[C@@H](O)[C@@H](O)[C@H]1CO
    NeuAc,O1[C@@H]([C@H](O)[C@H](O)CO)[C@H](NC(C)=O)[C@@H](O)CC1(O)C(=O)O
    Ne{u,

Verbosity
-----

As we see form the code-cell above, the program might print some messages. This can be suppressed on different levels using the built-in python module logging.

.. code:: python

    convert(glycan_list=["Gal", "Man", "GalNAc", "NeuAc", "Ne{u"], verbose=logging.WARNING)

Output:

.. code:: console

    ERROR:root:A parsing error occurred with Ne{u: <class 'glyles.glycans.utils.ParseError'>
    Error message: Glycan cannot be parsed:
    line 1:2 token recognition error at: 'e{'
    [('Gal', 'O1C(O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('Man', 'O1C(O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1CO'),
     ('GalNAc', 'O1C(O)[C@H](NC(C)=O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('NeuAc', 'O1[C@@H]([C@H](O)[C@H](O)CO)[C@H](NC(C)=O)[C@@H](O)CC1(O)C(=O)O'),
     ('Ne{u', '')]

As pythons built-in logging library does not support for turning off all loggings, we implemented this on our own for the input verbose=None.

.. code:: python

    convert(glycan_list=["Gal", "Man", "GalNAc", "NeuAc", "Ne{u"], verbose=None)

Output:

.. code:: python

    [('Gal', 'O1C(O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('Man', 'O1C(O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1CO'),
     ('GalNAc', 'O1C(O)[C@H](NC(C)=O)[C@@H](O)[C@@H](O)[C@H]1CO'),
     ('NeuAc', 'O1[C@@H]([C@H](O)[C@H](O)CO)[C@H](NC(C)=O)[C@@H](O)CC1(O)C(=O)O'),
     ('Ne{u', '')]