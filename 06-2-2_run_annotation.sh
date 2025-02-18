#!/bin/bash

# load the libraries
source /exchange/healthds/singularity_functions
module load python-3.9.10/py-requests/2.31.0
module load python-3.9.10/py-pandas/1.5.3

# run the VEP annotator
# replace 'proxy_snps' with the desired output prefix for your VEP annotated file
python vep_annotator.py snps_for_annotation_163.txt --output-prefix proxy_snps
