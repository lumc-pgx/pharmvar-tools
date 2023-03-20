# Introduction

Set of tools to work with pharmvar data.

# Installation (Linux)

```sh
python3 -m venv
source venv/bin/activate
pip install -r requirements.txt
```

# Installation (Windows)

1. Install Windows subsystem for Linux (WSL).

```sh
wsl --install
```

2. Navigate to the repository and install the Linux/Python dependencies.

```sh
cd /mnt/c/[path]/pharmvar-tools

sudo apt update
sudo apt upgrade
sudo apt install python3.10-venv python-is-python3 python3-pip

sudo python -m venv venv
source ./venv/bin/activate
python -m pip install -r requirements.txt
```

# Data requirements

For each gene of interest the sequence fasta of the NC and the NG should be
downloaded (from the NCBI) into the `data/` directory, e.g.,
`data/NC_000022.11.fasta` and `data/NG_008376.4.fasta` for *CYP2D6*.

In addition, for the check functionality, also the complete download zipfile
from pharmvar is required. It needs to be unzipped in the data directory, e.g., 
the fasta file for *CYP2D6* should be located in
`data/pharmvar-5.2.19/CYP2D6/CYP2D6.haplotypes.fasta`.

# Usage

## Pharmvar data consistency check

`python check.py --gene CYP2D6 --all`


## Calculating relations between alleles and variants based on variant algebra

`python compare.py --gene CYP2D6 > data/pharmvar_5.2.19_CYP2D6_relations.txt`


## Visualizing relations

Requires graphviz (`sudo apt install graphviz`).

`python to_dot.py --gene CYP2D6 < data/pharmvar_5.2.19_CYP2D6_relations.txt | fdp -Tpdf > /tmp/CYP2D6.pdf`
