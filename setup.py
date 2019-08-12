import os
import re
from setuptools import setup

NAME = 'isb_firm'
PACKAGES = ['firm']
DESCRIPTION = 'Framework for inference  of regulation by miRNAs'
LICENSE = 'LGPL V3'
URI = 'https://github.com/baliga-lab/firm2'
AUTHOR = 'Baliga Lab, Institute for Systems Biology'
VERSION = '1.0.0'

KEYWORDS = ['miRvestigator', 'firm']

# See trove classifiers
# https://testpypi.python.org/pypi?%3Aaction=list_classifiers

CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules"
    ]

PACKAGE_DATA = {
    'firm': ['common/*', 'FreqFiles/*',
             'TargetPredictionDatabases/PITA/*',
             'TargetPredictionDatabases/TargetScan/*']
}

if __name__ == '__main__':
    setup(name=NAME, description=DESCRIPTION,
          license=LICENSE,
          url=URI,
          version=VERSION,
          author=AUTHOR,
          author_email='wwu@systemsbiology.net',
          maintainer=AUTHOR,
          maintainer_email='wwu@systemsbiology.net',
          keywords=KEYWORDS,
          packages=PACKAGES,
          zip_safe=False,
          classifiers=CLASSIFIERS,
          include_package_data=True, package_data=PACKAGE_DATA,
          scripts=['bin/firm', 'bin/firm-convertminer'])
