Informatics Matters Standardize_smiles Utility
==============================================

The ``Standardize_smiles`` module is a utility provided by `Informatics Matters`_' computational pipelines.
It contains a python script to take an original SMILES (*Simple Molecule Input Line rESpresentation*) and:
translate it to an Isomeric or Nonisomeric cononical SMILES representation in a consistent way.

Preparation
-----------
Note that although this module is available on Pypi, it is not possible to simply install it and run it.
The module is strongly dependent on `RDKit`_.
In most cases, the Informatics Matters cartridge will be sufficient. More details can be found here: `IM-RDKit`_

Testing
-------

To check installation is successful, run::

    python -m unittest tests.tests

.. _Informatics Matters: http://www.informaticsmatters.com
.. _RDKit: https://www.rdkit.org/docs/index.html
.. _IM-RDKit: https://github.com/InformaticsMatters/rdkit_cartridge
