#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Setup module for the SMILES Utilities.
#
# March 2021

import os
from setuptools import setup, find_packages


def get_long_description():
    return open('README.rst').read()


setup(

    name='im-standardize-molecule',
    version=os.environ.get('GITHUB_REF_SLUG', '0.0.1'),
    author='Duncan Peacock',
    author_email='dpeacock@informaticsmatters.com',
    url='https://github.com/InformaticsMatters/standardize-molecule',
    license='GPLv3 License',
    description='Utilities for Molecular Science',
    long_description=get_long_description(),
    keywords='smiles, rdkit',
    platforms=['any'],

    # Our modules to package
    packages=find_packages(exclude=['*.test', '*.test.*', 'test.*', 'test']),
    py_modules=['standardize_molecule'],

    # Supported Python versions
    python_requires='>=3, <4',

    # Project classification:
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],

    zip_safe=False,

)
