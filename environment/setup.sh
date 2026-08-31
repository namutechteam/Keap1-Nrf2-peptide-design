#!/bin/bash
# Create the two conda environments used by this repository.
#
#   bash setup.sh
#
# Two are needed because chai_lab pins numpy 1.x and haddock3 pins numpy 2.x.
# Stages 01-04 run in `keap1nrf2`, stage 05 in `keap1nrf2-haddock`; they exchange
# PDB files on disk, so they never need to be importable together.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONNOUSERSITE=1

conda env create -f "$HERE/keap1nrf2.yml"
conda env create -f "$HERE/keap1nrf2-haddock.yml"

# plip runs `pip install requests` from its setup.py, which pip's build isolation
# blocks; requests is already installed from the yml, so turn isolation off.
conda run -n keap1nrf2 --no-capture-output \
    pip install --no-build-isolation "plip==2.4.0"

echo "done"
