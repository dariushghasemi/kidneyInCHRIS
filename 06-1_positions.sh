#!/bin/sh

LD=~/projects/gwas/pairways_LD/LD_range/LD

head -1 $LD/chr_1_10670853_LD.txt > test_convert.txt

ls -v $LD | while read line; do
	chr=$(echo $line | cut -d_ -f2)
	snp=$(echo $line | cut -d_ -f3) 
	pos=$(echo chr${chr}:${snp})
	cat $LD/${line}  | awk "{ if( ((\$2==$snp) || (\$3==$snp)) && (\$5 > 0.8) ) {print } }" >> test_convert.txt
done

cat test_convert.txt | sed -re 's/chr//' | awk ' NR > 1 {print $1"\t"$2"\n"$1"\t"$3}' | sort -g -k 1,2 | uniq > test_convert_prepared.txt
