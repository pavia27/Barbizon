#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

chmod +x barbizon.py

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

## adding conda channels
conda config --add channels defaults 2> /dev/null
conda config --add channels bioconda 2> /dev/null
conda config --add channels conda-forge 2> /dev/null
conda config --add channels au-eoed 2> /dev/null

## creating GToTree environment and installing dependencies
conda create -n barbizon -c r r-argparse pandas numpy prodigal scikit-learn metabat2 --yes

## activating environment and adding in kaiju if it didn't work before
source activate barbizon
conda install -c bioconda cat --yes

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding path to executable script
export PATH=\"$(pwd):"'$PATH'\"" \

echo '#!/bin/sh'" \
export rscripts=\"$(pwd)\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# re-activating environment so variable and PATH changes take effect
source activate barbizon


printf "\n        ${GREEN}DONE!${NC}\n\n"
