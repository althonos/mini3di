# `mini3di` [![Stars](https://img.shields.io/github/stars/althonos/mini3di.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/mini3di/stargazers)

*A [NumPy](https://numpy.org/) port of the [`foldseek`](https://github.com/steineggerlab/foldseek) code for encoding structures to 3di.*
<!-- 
[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/mini3di/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/mini3di/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/mini3di?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/mini3di/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/mini3di.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/mini3di)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/mini3di?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/mini3di)
[![Wheel](https://img.shields.io/pypi/wheel/mini3di.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/mini3di/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/mini3di.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/mini3di/#files)
[![Python Implementations](https://img.shields.io/badge/impl-universal-success.svg?style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/mini3di/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/mini3di/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/mini3di/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/mini3di.svg?style=flat-square&maxAge=600)](https://github.com/althonos/mini3di/issues)
[![Docs](https://img.shields.io/readthedocs/mini3di/latest?style=flat-square&maxAge=600)](https://mini3di.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/mini3di/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fmini3di)](https://pepy.tech/project/mini3di) -->

## 🗺️ Overview

`mini3di` is a pure-Python package to encode 3D structures of proteins into 
the 3di alphabet, using the trained weights from the `foldseek` VQ-VAE model.

This library only depends on NumPy and is available for all modern Python
versions (3.7+). 

### 📋 Features


### 🧊 Vectorization


## 🔧 Installing

For now, the package can be installed from the GitHub repository with `pip`:
```console
$ pip install git+https://github.com/althonos/mini3di
```

<!-- Install the `mini3di` package directly from [PyPi](https://pypi.org/project/mini3di)
which hosts universal wheels that can be installed with `pip`:
```console
$ pip install mini3di
``` -->

<!-- Otherwise, `mini3di` is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda mini3di
``` -->

<!-- ## 📖 Documentation

A complete [API reference](https://mini3di.readthedocs.io/en/stable/api.html)
can be found in the [online documentation](https://mini3di.readthedocs.io/),
or directly from the command line using
[`pydoc`](https://docs.python.org/3/library/pydoc.html):
```console
$ pydoc mini3di
``` -->

## 💡 Example

`mini3di` provides a single `Encoder` class, which expects the 3D coordinates
of the **Cα**, **Cβ**, **N** and **C** atoms from each peptide residue. For 
residues without **Cβ** (Gly), simply write the coordinates as `math.nan`. 
Call the `encode_atoms` method to get a sequence of 3di states:
```python
from math import nan
import mini3di

encoder = mini3di.Encoder()
states = encoder.encode_atoms(
    ca=[[32.9, 51.9, 28.8], [35.0, 51.9, 26.6], ...],
    cb=[[ nan,  nan,  nan], [35.3, 53.3, 26.4], ...],
    n=[ [32.1, 51.2, 29.8], [35.3, 51.5, 28.1], ...],
    c=[ [34.4, 51.7, 29.1], [36.1, 51.1, 25.8], ...],
)
```

The states returned as output will be a NumPy array of state indices. To turn
it into a sequence, use the `build_sequence` method of the encoder:
```python
sequence = encoder.build_sequence(states)
print(sequence)
```

The encoder can work directly with Biopython objects, if Biopython is available.
A helper method `encode_chain` to extract the atom coordinates from
a [`Bio.PDB.Chain`](https://biopython.org/docs/latest/api/Bio.PDB.Chain.html) 
and encoding them directly. For instance, to encode all the chains from a 
[PDB file](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)):
```python
import pathlib

import mini3di
from Bio.PDB import PDBParser

encoder = mini3di.Encoder()
parser = PDBParser(QUIET=True)
struct = parser.get_structure("8crb", pathlib.Path("tests", "data", "8crb.pdb"))

for chain in struct.get_chains():
    states = encoder.encode_chain(chain)
    sequence = encoder.build_sequence(states)
    print(chain.get_id(), sequence)
```

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/mini3di/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/mini3di/blob/main/CONTRIBUTING.md)
for more details.

## ⚖️ License

This library is provided under the 
[GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/). 
It includes some code ported from `foldseek`, which is licensed under the 
[GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/)
as well.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original `foldseek` authors](https://github.com/steineggerlab). 
It was developed by [Martin Larralde](https://github.com/althonos/) during his 
PhD project at the [European Molecular Biology Laboratory](https://www.embl.de/) 
in the [Zeller team](https://github.com/zellerlab).*
