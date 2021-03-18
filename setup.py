#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Setup module for the Python-based Jenkins Utilities.
#
# March 2021

import platform
from setuptools import setup, find_packages


def get_long_description():
    return open('README.rst').read()


setup(

    name='im-standardize-smiles',
    version='0.0.1',
    author='Alan Christie',
    author_email='achristie@informaticsmatters.com',
    url='https://github.com/InformaticsMatters/jenkins-utils',
    license='Copyright (C) 2021 Informatics Matters Ltd. All rights reserved.',
    description='Utilities for Informatics Matters CI/CD configuration',
    long_description=get_long_description(),
    keywords='smiles, rdkit',
    platforms=['any'],

    # Our modules to package
    packages=find_packages(exclude=['*.test', '*.test.*', 'test.*', 'test']),
    py_modules=['standardize_smiles'],

    # Supported Python versions
    python_requires='>=3, <4',

    # Project classification:
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Topic :: Software Development :: Build Tools',
    ],

    zip_safe=False,

)
