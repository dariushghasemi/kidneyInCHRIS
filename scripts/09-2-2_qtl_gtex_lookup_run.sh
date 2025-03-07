#!/bin/bash

#SBATCH --job-name gtex10
#SBATCH --output %j_query_gtex10.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --time 30-00:00:00

# load the libraries
source /exchange/healthds/singularity_functions
module load python-3.9.10/py-requests/2.31.0
module load python-3.9.10/py-pandas/1.5.3
module load python-3.9.10/py-tqdm/4.66.1

# run query in GTEx rest API
python query_gtex.py \
    --variants_file mediation_table2.variants \
    --output_prefix "02-Mar-25" \
    --dataset 'gtex_v10' \
    --pvalue 0.1 \
    --threads 3
