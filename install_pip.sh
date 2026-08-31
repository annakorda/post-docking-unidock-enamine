#!/bin/bash
# Author: Anna Korda
# Run this with the luna-chem conda env already activated (see environment.yml).
# Order matters here -- verified from a genuinely clean pip cache (no
# pre-built wheels to fall back on), not just re-tested with a cache that
# masked the real problem the first time around:
#
#   1. luna's own declared runtime deps first (colorlog, xopen, mmh3,
#      PDBeCif) -- matches `pip show luna`'s Requires: line.
#   2. dimorphite-dl, meeko, gemmi, loguru -- installing dimorphite-dl
#      here (not after) is what actually supplies rdkit, which luna's
#      setup.py imports at build time but never declares as a dependency.
#   3. luna itself, with --no-build-isolation, now that everything its
#      setup.py needs to import (networkx, rdkit, ...) is already present
#      in this env instead of an empty pip build sandbox.
set -e

pip install colorlog xopen "mmh3==2.5.1" "PDBeCif==1.5"
pip install "meeko==0.7.1" gemmi dimorphite-dl loguru
pip install --no-build-isolation "luna==0.14.0"

python -c "
import luna, meeko, dimorphite_dl, loguru, rdkit, networkx
from pymol import cmd
from openbabel import openbabel
import Bio
print('luna+chemtools env: all imports OK')
"
