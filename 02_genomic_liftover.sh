#!/bin/sh

BASE=/home/dghasemisemeskandeh/projects/gwas/Data/CHRIScohort
CHAIN=/home/demmert/work/Functionalization/work/hg19ToHg38.over.chain.gz
META=$BASE/meta.results_corrected.with.MetaSubtract.txt.gz
OUT=$BASE/02_liftover_ckdgen/output

#--------------#
# First we take the genomic cooardinates
# from subtracted CKDGen summary stats

# take CHR and POS clumns
zcat  ${META} | cut -f1,2 > ${OUT}/POS.txt

# modify the chromosomes number
cat ${OUT}/Bed.txt | sed -e 's/^/chr/g ' -e 's/chrChr/CHR/g' | tail -n+2 > ${OUT}/positions.bed

#------------------------------------------#
#    lift to hg38 via CrossMap v0.5.3      #
#------------------------------------------#

# perform lift-over
CrossMap.py bed ${CHAIN}
            ${OUT}/positions.bed
	    ${OUT}/positionsHg38.bed

#--------------#
# Joining  snps lifted positions to coordinates in build 37 
cat ${OUT}/positionsHg38.bed | cut -f1,2,5,6 | sed -e 's/chr*//g' >> ${OUT}/MetaBED38.txt

#--------------#
# run on the servers
# sarrayscript  -c 8 --mem-per-cpu=8GB -J 02_liftover.sh

 