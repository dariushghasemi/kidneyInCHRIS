#!/bin/bash

BASE=/home/dghasemisemeskandeh/projects/gwas
PHENO=$BASE/Data/phenotypes/phenodf4_NEW_scheme_W.ped
OUT=$BASE/Output/ReformedScheme
KIN=/scratch/ekoenig/CHRIS_CORRECTED/10K/Genotyped/20210409/CHRIS.10K.genotypes.pruned.MAF0.01.kinf
DATA=/scratch/ekoenig/CHRIS_CORRECTED/10K/Imputed/TOPMedR2/20210409
BIN=/home/cfuchsberger/bin/epacts326.topmed/bin/


for t in SerumCreatinine.Stdw.Res	eGFRw.Res	eGFRw.log.Res
do
  for i in `seq 1 22`
  do
     echo "${BIN}epacts single -vcf $DATA/chr$i.vcf.gz \
          -ped ${PHENO} --pheno $t \
          -out ${OUT}/$t.chr$i \
          -test q.emmax \
          -kinf ${KIN} --chr $i \
          -field DS \
          --run 24 --mosix-nodes \"\""
  done
done
																											    