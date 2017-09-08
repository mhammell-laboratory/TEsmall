#!/usr/bin/env python

from setuptools import setup, find_packages

execfile('tesmall/version.py')

with open("README.md") as f:
    readme = f.read()

with open('requirements.txt') as f:
    install_requires = f.read().strip().splitlines()

setup(name="TEsmall",
      version=__version__,
      description="A pipeline for profiling TE-associated small RNAs",
      long_description=readme,
      classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      keywords="TE transposable element small RNA",
      url="https://github.com/wwliao/tesmall",
      author="Wen-Wei Liao",
      author_email="gattacaliao@gmail.com",
      license="GPLv3",
      packages=find_packages(),
      install_requires=install_requires,
      entry_points={
        'console_scripts': ["tesmall=tesmall.command_line:main"],
      },
      include_package_data=True,
      zip_safe=False)
