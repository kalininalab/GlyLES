Installation
============


PIP installation
----------------

So far, this package can only be downloaded from the python package index. So the installation with ``pip`` is very easy.
Just type

.. code:: console

    pip install glyles


and you're ready to use it as described below. Use

.. code:: console

    pip install --upgrade glyles

to upgrade the glyles package to the most recent version.

Testing installation
--------------------

To test the installation, just run

.. code:: python

    from glyles import convert

    convert("Gal")

and it should output a SMILES string for galactose.