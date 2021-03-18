#!/usr/bin/env python
# coding=utf-8

"""standardize_molecule.py

Utility to standardize an input smiles into the standard representation
used in the fragment network.

The module provides methods for rendering a vendor (original) SMILES
representation into a standard smiles format in isomeric or nonisomeric form
The nonisomeric form is used for searching within the fragment network.

Dependencies
------------
Note that this is heavily dependent on RDKit.

Duncan Peacock
March 2021
"""

import logging

# RDKit
# pylint: disable=E0401
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# Our logger
logger = logging.getLogger(__name__)


def get_biggest_component(mol):
    """Get largest component molecule.
    """
    frags = Chem.GetMolFrags(mol, asMols=True)

    if len(frags) == 1:
        return mol

    i = 0
    biggest_count = 0
    for frag in frags:
        hac = frag.GetNumHeavyAtoms()
        if hac > biggest_count:
            biggest_count = hac
            biggest_mol = frag
        i += 1

    return biggest_mol


def remove_isotopes(mol):
    """Remove isotopes from molecule
    """
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)


def standardize_format(mol):
    """Clean up molecule and return in standardized format
    """
    mol = rdMolStandardize.Cleanup(mol)
    mol = get_biggest_component(mol)

    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    remove_isotopes(mol)
    return mol


def standardize_molecule(mol, osmiles = 'NotProvided'):
    """Standardise a RDKit mol object. A SMILES string that represents the
    molecule can also be passed in to allow better logging of errors.

    :param mol: flag indicating whether this is a molecule or not.
    :param osmiles: The original (non-standard) SMILES

    :return: standardized form of molecule
             if error the standard form will be returned as None.
     """

    std = None
    try:
        std = standardize_format(mol)
    except Exception as ex:
        logger.warning('standardize(%s) exception: "%s"',
                       osmiles, ex.message)

    if not std:
        logger.error('Got nothing from standardize(%s).'
                     ' Skipping this compound', osmiles)
    return std


def get_smiles(isomeric_smiles, std_molecule, osmiles):
    """Return Smiles of desired format given a standardised molecule
    """

    smiles = None
    try:
        smiles = Chem.MolToSmiles(std_molecule, isomeric_smiles, canonical=True)
    except Exception as ex:
        logger.warning('MolToSmiles(%s, noniso) exception: "%s"',
                       osmiles, ex.message)
    if not smiles:
        logger.error('Got nothing from MolToSmiles(%s, noniso).'
                     ' Skipping this compound', osmiles)
    return smiles


def process_osmiles(osmiles, isomeric_smiles):
    """Takes original smiles and processes to return a smiles of the
    desired format (isomeric or nonisomeric) with the standardized molecule
    """
    mol = None
    try:
        mol = Chem.MolFromSmiles(osmiles)
    except Exception as ex:
        logger.warning('MolFromSmiles(%s) exception: "%s"',
                       osmiles, ex.message)

    smiles = None
    std_molecule = None
    if mol:
        std_molecule = standardize_molecule(mol, osmiles)
        if std_molecule:
            smiles = get_smiles(isomeric_smiles, std_molecule, osmiles)
    else:
        logger.error('Got nothing from MolFromSmiles(%s).'
                     ' Skipping this compound', osmiles)

    return smiles, std_molecule


def standardize_to_iso_smiles(osmiles):
    """Given a vendor (original) SMILES this method standardises
    it into a canonical form and returns a tuple that contains
    the isomeric form and the standard form.

    :param osmiles: The original (non-standard) SMILES

    :returns: Standardized isometric SMILES representation
              Standardized molecule.
    On error, the standard form will be returned as None.
    """
    if not osmiles:
        return None

    return process_osmiles(osmiles, True)


def standardize_to_noniso_smiles(osmiles):
    """Given a vendor (original) SMILES this method standardises
    it into a canonical form and returns a tuple that contains
    the non-isomeric form and the standard form.

    :param osmiles: The original (non-standard) SMILES

    :returns: Standardized nonisomol SMILES representation
              Standardized molecule.
    On error, the standard form will be returned as None.
    """
    if not osmiles:
        return None

    return process_osmiles(osmiles, False)
