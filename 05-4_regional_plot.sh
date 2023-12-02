#!/bin/sh

BASE=/home/dghasemisemeskandeh/projects/gwas
SNP=$BASE/Data/chrPos_LD.txt
GWAS=$BASE/Output/ReformedScheme/eGFRw.log.Res.txt.gz
LD=$BASE/04_pairways_LD/LD_range/LD_prepared
OUT=$BASE/05_regional_association/LZ_plots/filtered_eGFRw.log.Res_LZ_Dprime/

#---------------------------------------------------#
# The script generates regional association plots   #
# using stand alone LocusZoom based on pre-computed #
# local LD of the variants via emeraLD.             #
#---------------------------------------------------#

# Read positions of the leading variant at 
# each of the 147 kidney loci in CKDGen. 

for pos in $(cat $LD)
do
  chr=$(echo $pos | cut -d: -f1)
  fn=$(echo $pos | sed -e 's/:/_/')
  echo locuszoom --epacts ${GWAS}   \
                 --markercol MARKER_ID    \
		 --epacts-beg-col BEG    \
		 --epacts-end-col END   \
		 --epacts-chr-col CHROM    \
		 --chr ${chr}    \
		 --flank 500kb     \
		 --refsnp ${pos}    \
		 --ld ${LD}/chr_${fn}_LD_prepared.txt    \
		 --ld-measure dprime   \
		 --build hg38    \
		 --plotonly    \
		 --prefix '"filtered_eGFRw.log.Res_Dprime_211112_"'
done