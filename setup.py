#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import sys

packages = ['genomfart',
            'genomfart.parsers',
            'genomfart.data',
            'genomfart.utils',
            'genomfart.test',
            'genomfart.test.parsers',
            'genomfart.test.utils'
           ]
setup(
    name='genomfart',
    version='0.3',
    author='Eli Rodgers-Melnick',
    author_email='er432@cornell.edu',
    description='A Genomics package for Python',
    url='https://github.com/er432/genomfart',
    platforms=['Linux','Mac OSX', 'Windows', 'Unix'],
    keywords=['Genomics','parsers'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis'
    ],
    packages=packages,
    package_data={'genomfart':['data/*.gff']},
    license='BSD'
    )
