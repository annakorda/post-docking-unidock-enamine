#!/bin/bash
# Author: Anna Korda
# Run this with the post-dock conda env already activated (see environment.yml).
#
# --no-deps is load-bearing here, found the hard way: a plain
# `pip install unigbsa==0.1.7` pulls in "openbabel-wheel" as a transitive dep
# of acpype -- a second, pip-built openbabel that silently shadows the
# conda-installed one already in environment.yml and is ABI-incompatible with
# it (obabel then fails with "undefined symbol:
# _ZN9OpenBabel8OBPlugin7DisplayERSsPKcS3_", which breaks every ligand's
# sdf->mol conversion, i.e. every single compound). Every other real
# dependency unigbsa needs is already satisfied by environment.yml (acpype,
# gmx-mmpbsa, pandas, scipy, mpi4py) except lickit, which has no conda-forge
# package and is pure-python (no compiled extension, so no ABI risk) --
# installed explicitly instead.
set -e

pip install --no-deps "unigbsa==0.1.7" "lickit==0.1.0"

unigbsa-pipeline --help > /dev/null
gmx --version > /dev/null
pdb4amber --help > /dev/null
obabel -V > /dev/null
python -c "from rdkit import Chem; import numpy, matplotlib"
echo "post-dock env: unigbsa-pipeline, gmx, pdb4amber, obabel, rdkit, numpy, matplotlib all OK"
