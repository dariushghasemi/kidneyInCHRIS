<img src="https://github.com/DariushG3/kidneyInCHRIS/blob/main/Untitled.jpg" width=250 align="right">

## KidneyInCHRIS part of TrainCKDis Horizon 2020 EU Project
Genetic determinants of Kidney function in The **C**orporative **H**ealth **R**esearch **I**n **S**outh Tyrol (CHRIS) study

Several issues have been created to track the progress in the project. Here is the list of the issues in chronological order:

* [1- Picking the discovery Loci](#1-picking-the-discovery-loci)
* [2- Subtracting CHRIS from CKDGen](#2-subtracting-chris-from-ckdgen)
* [3- Phenotype creation](#3-phenotype-creation)
* [4- Conducting GWAS on renal function](#4-conducting-gwas-on-renal-function)
* [5- Comparing replicated variants effect magnitude with discovery studies](#)
* [6- Comparing the 147 CKDGen Loci in CHRIS](#)
* [7- Excluding batch effects from the health traits using transformation](#)
* [8- Extracting dosage from VCF](#8-extracting-dosage-from-vcf)
* [9- Phenome-wide mediation analysis](#9-phenome-wide-mediation-analysis)
* [10- Kidney interaction with Thyroid](#10-kidney-interaction-with-thyroid)
* [11- Sensitivity Analyses](#11-sensitivity-analyses)
* [12- Interrogating SNPs in Phenoscanner](#12-interrogating-snps-in-phenoscanner)
* [13- Publication and manuscript revision](#13-publication-and-manuscript-revision)
___________________________________________________________________________________________________________


## 1. Picking the discovery Loci
- We took 147 eGFR-associated loci identified by the CKDGen Consortium European ancestry GWAS meta-analysis in CKDGen (from Wuttke et al., Nat Genet 2019). One index SNP was available per locus in CKDGen. We then picked the proxy SNPs in the 1M bp region around the CKDGen index SNP with which were in strong LD (Rsq > 0.8) in the CHRIS study. We explain the replication process through GWAS in points 4 & 5.

## 2. Subtracting CHRIS from CKDGen
- We were going to indipendentizing CKDGen results by subtracting CHRIS 5K from CKDGen Meta-GWAS on European ancestry.

## 3. Phenotype creation

## 4. Conducting GWAS on renal function
- We mapped the genotypes of CHRIS 10K participants to our interested phenotype, __ln(eGFR)__, by conducting genome-wide association study (GWAS) analysis.

- To draw regional association plot for the lead variant of each replicated loci, we needed to update the LocusZoom SQLite database to show the recombination rate and it hotspots in the plot for SNPs position in build 38 (GRCh38). See more info in [issue 10](https://github.com/DariushG3/kidneyInCHRIS/issues).

## 5. Comparing replicated variants effect magnitude with discovery studies

## 6. Comparing the 147 CKDGen Loci in CHRIS

## 7. Excluding batch effects from the health traits using transformation

## 8. Extracting dosage from VCF
- To do association tests analysis, we extractied the dosage level of the 163 replicated SNPs from imputed VCF files. To do so, we used bcftools:
```bash

```

## 9. Phenome-wide mediation analysis

## 10. Kidney interaction with Thyroid




## 11. Sensitivity Analyses
### Sensitivity of SNPs' dosage to having missing in FT3/4
- Here we test if there is any difference between the dosage level for the replicated SNPs for those cases having missing in FT3 or FT4 due to the artifacts resulted from using a two-sided Wilcoxon test (alpha = 0.05):
```R
#Sensitivity test for missing in FT3:
wilcox.test(SNP ~ FT3_indicator, data = vcfReg)

#Sensitivity test for missing in FT4:
wilcox.test(SNP ~ FT4_indicator, data = vcfReg)
```
### Sensitivity of the FT3/4-adjusted association of SNPs on log(eGFR) to Municipality
- Here we test how the association of SNPs on log(eGFR) changes when we adjust the LR model both for FT3/4 and municipality. This variable indicates the unmeasured differences in allele frequency of the SNPs between different villages in the region where CHRIS study was conducted. Truly, this idea is coming from where we thought adding municipality could affect the SNP-Kidney adjusted association FT3 or FT4, and in this way we might capture the batch effects in order to have an unbiased idea of the interaction between the kidney and thyroid traits.
- Thus, we tested the association of the SNP and kidney adjusted for FT3/4 once without adjusting for municipality and once with adjusting for it as below:
```R
#Sensitivity test for FT3:
log(eGFR) ~ SNP +FT3 + Sex + Age vs. log(eGFR) ~ SNP +FT3 + Sex + Age + Municipality

#Sensitivity test for FT4:
log(eGFR) ~ SNP +FT4 + Sex + Age vs. log(eGFR) ~ SNP +FT4 + Sex + Age + Municipality
```

## 12. Interrogating 163 replicated SNPs in Phenoscanner




## 13. Publication and manuscript revision



- List of in-LD variants created using `06-1-2_ld_analysis_result.R` (25-Nov-23).

- Extraction of the dosage of in-LD variants was done in five minutes on the Eurac servers using `06-1-3_ld_variants_extraction.sh` (Mon, 19:05, 27-Nov-23).

- Plotted the variance of i-LD variants explained by each PCs `06-1-4_ld_variants_pca.R`. Still need to plot the number of PCs explaining different portion of the total variance of the data (Wed, 18:34, 29-Nov-23).

- Histogram and cumulative distribution of number of PCs explaining different proportions of the total variance of the in-LD variants was genrated and the report was written and sent to the PI (Thu, 20:04, 30-Nov-23).

- The scripts written to be used for running GWAS on ln(eGFRcrea) and doing replication study in the CHRIS study were added to Github (Sat, 23:55, 02-Dec-23).

- The variance explained by the 147 GWAS hits were quantified and compared in both CHRIS and CKDGen. The quantities were illustrated in proper charts via `07-3_variance_explained.R` (Sun, 19:35, 03-Dec-23). 

- Unnecessary scripts like `CHRIS.R` were removed from git repository and stored in local computer (Sun, 23:35, 03-Dec-23). 

- Scripts order were corrected again. The script to generate the qq-plot have been stored to be sent to the Reviewer (Mon, 13:40, 04-Dec-23).

- The amount the variance explained (97.5%) by the first 147 principal components were plotted as supplementary figure S1 using `06-1-4_ld_variants_pca.R` (Tue, 05-Dec-23).

- Replication workflow deposited in `06-2_extraction_loci_all.R` and `06-3_extraction_loci_replicated.R` was rephrase in tidyverse style for its reproducibility and DRY manner (Thu, 04:30, 07-Dec-23).

- Trans-ethnic CKDGen meta-GWAS analysis after excluding CHRIS 5K was added to the supple. Table 2. (Thu, 18:25, 07-Dec-23).

- The plot showing variants' effect in CKDGen vs. CHRIS (Fig 3B) made via `07-2_chris_vs_ckdgen_plots.R` was merged with the newly drawn plot showing variance explained by each of the lead variants at 147 kidney loci (Fig. 3A) made here `07-3_variance_explained.R`; now it is Fig. 3 of the revised paper (Sat, 03:34, 09-Dec-23).

- The 11 hits in paper Table 2 resulted from subtracting CHRIS 5K from CKDGen trans-ethnic meta-GWAS was saved below - rs-numbers were extracted from Table 2 of the revised paper (Sat, 04:04, 09-Dec-23).
```bash
CKDGen_TA="~/projects/gwas/01_subtract_ckdgen/output/meta.results_corrected.with.MetaSubtract.txt.gz"
OUTPUT="~/projects/gwas/06_replication_analysis/output/09-Dec-23_chris_11_rep_loci_in_subtracted_ckdgen_trans-ethnic.txt"

zgrep -E 'rs74748843|rs807624$|rs28817415|rs10062079|rs3812036|rs57514204|rs819196|rs2039424|rs7113042|rs59646751|rs77924615' $CKDGen_TA > $OUTPUT
```

- I got sick, probably covid-19 infection and couldn't work on the shroom3 draft paper. (Sat, 20:55, 09-Dec-23) 

- Just recovered and after submitting paper revision tomorrow noon, I will work on the 2nd paper to finish the discussion (Sun, 23:22, 10-Dec-23).

- Fig. 3A was modified and Fig.3 was changed accordingly(Mon, 16:30, 11-Dec-23).

- We submitted the revised manuscript to the journal (Mon, 23:59, 11-Dec-23).


Dariush

