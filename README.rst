Informatics Matters Standardize SMILES Utility
==============================================
The ``im-standardize_smiles`` module is a utility provided by `Informatics Matters`_.
It contains a Python module to take an original SMILES (*Simple Molecule Input Line rESpresentation*) and
translate it to an Isomeric or Nonisomeric cononical SMILES representation in a consistent way.

Preparation
-----------
Note that although this module is available on PyPI, it is not possible to
simply install it and run it. The module is strongly dependent on `RDKit`_.
In most cases, the Informatics Matters cartridge will be sufficient.
More details can be found here: `IM-RDKit`_

Testing
-------
To check installation is successful, from within an RDKit environment, run::

    python -m unittest test.test

.. _Informatics Matters: http://www.informaticsmatters.com
.. _RDKit: https://www.rdkit.org/docs/index.html
.. _IM-RDKit: https://github.com/InformaticsMatters/rdkit_cartridge
