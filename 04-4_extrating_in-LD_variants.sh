#!/usr/bin/bash

# base directory
BASE="/home/dghasemisemeskandeh/projects/kidneyInCHRIS"

# tab-separated list of 6337 in-LD variants (r2>0.8, build GRCh38)
variants_file="25-Nov-23_variants_in_ld.txt"

# TOPMedR2-imputed genotype data for CHRIS 10K (build 38)
genotype="/shared/statgen/CHRIS10K/Imputation/TOPMedR2"

# create the headers of the dosage file
echo -e "AID\tCHR\tPOS\tSNP_ID\tDosage" > $BASE/27-Nov-23_variants_in_ld_dosage.txt

# extracting the in-LD variants
for i in `seq 1 22`;
do
  bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%DS\n]' $genotype/chr$i.vcf.gz -R $variants_file >> $BASE/27-Nov-23_variants_in_ld_dosage.txt
done

# sbatch ./05-3_extrating_in-LD_variants.sh  -c 4 --mem-per-cpu=8GB -J "05-3.sh"

