#!/bin/bash

#SBATCH --job-name decode
#SBATCH --output %j_query_decode.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 30-00:00:00
#SBATCH --array=1-104

#-------------#
# Here we look-up of mediation variants in decode study 
# on proteomics (n=5141) to identify pQTLs.

# This script can be used with job schedulers like SLURM
# Example SLURM usage:
# sbatch --array=1-3 query_decode_run.sh

#-------------#
# Configuration
DATA_DIR="/exchange/healthds/public_data/sumstats/decode/largescaleplasma-2023/final_somascan_smp"
VARIANT_FILE="/scratch/dariush.ghasemi/projects/decode/mediation_table3.rsid"
OUTPUT_DIR="results"

#-------------#
# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# load the libraries
source /exchange/healthds/singularity_functions

#-------------#
# Determine the job array ID and total number of tasks
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    TASK_ID=$SLURM_ARRAY_TASK_ID
    TOTAL_TASKS=$SLURM_ARRAY_TASK_COUNT
else
    # For local testing
    TASK_ID=${1:-1}
    TOTAL_TASKS=${2:-10}
fi

#   zgrep -E 'rs28394165|rs10025351|rs28817415|rs4859682|rs13146355|rs3812036' "$path" >> ${OUTPUT_DIR}/results_10.tsv

#-------------#
# Read the list of desired variants into an array
VARIANTS=()
while IFS= read -r variant; do
  VARIANTS+=("$variant")
done < "$VARIANT_FILE"

#-------------#
# Create pattern for grep
echo "Creating grep pattern from $VARIANT_FILE"
GREP_PATTERN=$(paste -sd "|" "$VARIANT_FILE")

echo "Running job $TASK_ID of $TOTAL_TASKS"

#-------------#
# Get the list of all files
echo "Finding proteomics files in $DATA_DIR..."
mapfile -t all_files < <(find "$DATA_DIR" -name "Proteomics_SMP_PC0_*.txt.gz" | sort)
total_files=${#all_files[@]}
echo "Found $total_files total files"

#-------------#
# Calculate which files this task should process
files_per_task=$(( (total_files + TOTAL_TASKS - 1) / TOTAL_TASKS ))
start_idx=$(( (TASK_ID - 1) * files_per_task ))
end_idx=$(( start_idx + files_per_task - 1 ))

if [ $end_idx -ge $total_files ]; then
    end_idx=$(( total_files - 1 ))
fi

#-------------#
# output filename
output_file="$OUTPUT_DIR/results_${TASK_ID}.tsv"
echo "Results will be saved to: $output_file"

#-------------#
# Create a header for the output file
echo -e "Chrom\tPos\tName\trsids\teffectAllele\totherAllele\tBeta\tPval\tminus_log10_pval\tSE\tN\tImpMAF\tseqid\tgene_name\tprotein_name\tfilename" > "$output_file"

# show progress
echo "Task $TASK_ID processing files $start_idx to $end_idx (out of $total_files)"

#-------------#
# Process each file assigned to this task
for (( i=start_idx; i<=end_idx; i++ )); do
    file="${all_files[$i]}"
    filename=$(basename "$file")
    
    # Show progress every 10 files
    if (( (i - start_idx) % 10 == 0 )); then
        echo "Processing file $(( i - start_idx + 1 )) of $(( end_idx - start_idx + 1 ))"
    fi
  
  # Extract metadata from filename
  # Example: Proteomics_SMP_PC0_10000_28_CRYBB2_CRBB2_10032022.txt.gz
  if [[ $filename =~ Proteomics_SMP_PC0_([0-9]+_[0-9]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)_ ]]; then
    seqid="${BASH_REMATCH[1]}"
    gene_name="${BASH_REMATCH[2]}"
    protein_name="${BASH_REMATCH[3]}"
    
    # Process the file - skip header, grep for variants
    zgrep -E $GREP_PATTERN $file |\
    awk -v seqid="$seqid" -v gene="$gene_name" -v protein="$protein_name" -v fname="$filename" '
      BEGIN { OFS="\t" }
      NR > 1 { 
          print $0, seqid, gene, protein, fname; 
          }
      ' >> $output_file
    
    # Report how many matches were found in this file
    echo "  Found matches in $filename"
  else
    echo "Warning: Could not parse metadata from filename: $filename" >&2
  fi
done

#-------------#
# Report final counts
total_matches=$(wc -l < "$output_file")
total_matches=$((total_matches - 1))  # Subtract header line

echo "Task $TASK_ID completed - Found $total_matches total matches across all processed files"

#-------------#
# Flush any pending output
sleep 1

exit 0