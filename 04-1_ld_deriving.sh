#!/bin/sh


BASE=/home/dghasemisemeskandeh/projects/gwas
SNP=$BASE/Data/chrPos_LD.txt
DATA=/scratch/ekoenig/CHRIS_CORRECTED/10K/Imputed/TOPMedR2/20210409
tail=500000

for pos in $(cat $LD)
do
	chr=$(echo $pos | cut -d: -f1)
	fn=$(echo $pos | sed -e 's/:/_/')
	snp=$(echo $pos | cut -d: -f2)
	L="$(($snp-$tail))"
	U="$(($snp+$tail))" #| bc
	echo emeraLD -i  $DATA/chr$chr.vcf.gz --region chr$chr:$L-$U --dstats --stdout \> chr_${fn}_LD.txt
done


