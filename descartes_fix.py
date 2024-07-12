#!/bin/bash
conda env create -f environment.yml
source activate astropy-env
patch $(conda env list | grep astropy-env | awk '{print $2}')/lib/python3.9/site-packages/patch.py < descartes_fix/patch.py