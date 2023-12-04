#!/bin/sh

BASE=/home/dghasemisemeskandeh/projects/gwas
SNP=$BASE/Data/chrPos_LD.txt
OUTPUT=$BASE/pairways_LD/LD_range

for pos in $(cat $SNP)
do
	chr=$(echo $pos | cut -d: -f1)
	fn=$(echo $pos | sed -e 's/:/_/')
	
	echo    cat $OUTPUT/LD/chr_${fn}_LD.txt '|'  \
		sed -r "'s/chr//g'" '|'   \
		awk -F '"\\t"' $(echo '{t= $1; $1= t":"$3; $2= t":"$2; $3= $NF; $4=$5; print ;}')  OFS='"\\t"' '|'   \
		cut -f1-4 '|' sed -e "'1s/#CHR:POS2/snp1/'" -e "'1s/#CHR:POS1/snp2/'" -e "'1s/Dprime/dprime/'"  -e "'1s/Rsq/rsquare/'" \> $OUTPUT/LD_prepared/chr_${fn}_LD_prepared.txt
done
