#!/usr/bin/env python

from setuptools import setup, find_packages

exec(compile(open('TEsmall/version.py', "rb").read(), 'TEsmall/version.py', 'exec'))

with open("README.md") as f:
    readme = f.read()

with open('requirements.txt') as f:
    install_requires = f.read().strip().splitlines()

setup(name="TEsmall",
      version=__version__,
      description="A pipeline for profiling TE-associated small RNAs",
      long_description=readme,
      classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      keywords="TE transposable element small RNA",
      url="https://hammelllab.labsites.cshl.edu/software#TEsmall",
      author="Wen-Wei Liao, Kat O'Neill, Molly Hammell",
      author_email="koneill@cshl.edu",
      license="GPLv3",
      packages=find_packages(),
      install_requires=install_requires,
      entry_points={
        'console_scripts': ["TEsmall=TEsmall.command_line:main"],
      },
      include_package_data=True,
      zip_safe=False
)
