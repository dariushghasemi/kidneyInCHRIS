#!/usr/bin/bash


BASE=/home/dghasemisemeskandeh/projects/gwas
DATA=/scratch/ekoenig/CHRIS_CORRECTED/10K/Imputed/TOPMedR2/20210409
SNP=$BASE/Data/chrPos_LD.txt
OUT=$BASE/04_pairways_LD/LD_range/LD
tail=500000

for pos in $(cat $SNP)
do
  chr=$(echo $pos | cut -d: -f1)
  fn=$(echo  $pos | sed -e 's/:/_/')
  snp=$(echo $pos | cut -d: -f2)
  L="$(($snp-$tail))"
  U="$(($snp+$tail))"

  echo emeraLD -i  $DATA/chr$chr.vcf.gz \
               --region chr$chr:$L-$U  \ 
	       --dstats   \
	       --stdout \> $OUT/chr_${fn}_LD.txt

done


