Informatics Matters Standardize Molecule Utility
================================================
The ``im-standardize-molecule`` module is a utility provided by `Informatics Matters`_.

It contains a Python module to:
1. Take an original SMILES (*Simple Molecule Input Line rESpresentation*) and
convert it to an Isomeric or Nonisomeric canonical SMILES representation in a consistent way.
2. Convert an *rdkit.Chem.rdchem.Mol* object to a standard format.


Preparation
-----------
Note that although this module is available on PyPI, it is not possible to
simply install it and run it. The module is strongly dependent on `RDKit`_.
In most cases, one of the the Informatics Matters RDKit container images can be used.
More details can be found here: `IM-RDKit`_
RDKit can also be installed with `Conda`_.

Testing
-------
To check installation is successful, from within an RDKit environment, run::

    python -m unittest test.test

Usage
-----

**standardize_to_noniso_smiles** *(osmiles:str)*

**standardize_to_iso_smiles** *(osmiles:str)*

Takes a smiles and returns the converted smiles and standardized molecule as an *rdkit.Chem.rdchem.Mol*

Example::

    from standardize_molecule import standardize_to_noniso_smiles
    noniso = standardize_to_noniso_smiles('N[C@@H](C[C@@H]1CCNC1O)[C@@H](O)CO')
    [17:30:44] Initializing MetalDisconnector
    [17:30:44] Running MetalDisconnector
    [17:30:44] Initializing Normalizer
    [17:30:44] Running Normalizer
    [17:30:44] Running Uncharger
    print(noniso[0])
    NC(CC1CCNC1O)C(O)CO


**standardize_molecule** *(molecule: object, osmiles:str = 'NotProvided')*

Takes an rdkit.Chem.rdchem.Mol object and (optionally) a smiles string for logging and returns
a standardized molecule.

Example::

    from rdkit import Chem
    from standardize_molecule import standardize_molecule
    mol = Chem.MolFromSmiles('C1=CC=CC=C1')
    std = standardize_molecule(mol, 'C1=CC=CC=C1')
    [17:37:40] Initializing MetalDisconnector
    [17:37:40] Running MetalDisconnector
    [17:37:40] Initializing Normalizer
    [17:37:40] Running Normalizer
    [17:37:40] Running Uncharger
    print(std)
    <rdkit.Chem.rdchem.Mol object at 0x7f18ca2a93f0>
    print(noniso[1].GetNumHeavyAtoms())
    13


.. _Informatics Matters: http://www.informaticsmatters.com
.. _RDKit: https://www.rdkit.org/docs/index.html
.. _IM-RDKit: https://github.com/InformaticsMatters/docker-rdkit/blob/master/README.md
.. _Conda: https://www.rdkit.org/docs/Install.html#how-to-install-rdkit-with-conda