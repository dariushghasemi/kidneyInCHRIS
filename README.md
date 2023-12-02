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
- This analysis have been throughly explained earlier in this [repository](https://github.com/DariushG3/SNP-TSHcat_Interaction_Model)
- As we observed the significants interaction for two SNPs out of 3 in one replicated locus, we tried to elucidate this findings. By centralizing TSH (0-mean-centered) we realized adjustment in SNP regression coefficients but still interaction effects and the corresponding P-values remained stable.
- Another important thing we found was that the significant interaction only occurs when we exclude the cases with cancers/renal carcinoma/alternation of thyroid function due to pregnancy which were ~9,700 individuals. If we apply the previous model on the entire CHRIS 10K consisting of ~10,140 individuals the significance of interaction terms goes away and there would not be any significant interaction antymore.
```R
# Regular TSH levels
lm(eGFRw.log.Res ~ SNP * TSH + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
```

- We replaced log(eGFR) with eGFRw.log.Res in the above model and ran it again. The results stay consistent in a way that those two variants with significant inetraction with TSH remained significant. At the end, in this way we can avoid over-adjustment of the covariates: age and sex (03-Feb-23).

- We modified our exclusion criteria in the way that we omitted the observations with (09-Feb-2023):
1. missing values on both TSH and thyroid drugs (n = 4)
2. thyroid cancer (n = 16)
3. kidney cancer (n = 1)
4. goitre (n = 277)
5. operation on thyroid gland (n = 312)

- After excluding these individuals, exactly 9,730 people remined out of 10,146 for SNP-TSH interaction analysis (10-Feb-2023).
- Thus, the interaction models were executed as follows:
```R
# Interaction of SNP-TSH
lm(log(eGFR) ~ SNP * TSH     + PC1 + ... + PC10)
# Interaction of SNP-thyroid_disease
lm(log(eGFR) ~ SNP * TSH_cat + PC1 + ... + PC10)
```
- At the end, we found there was no significant interaction between kidney variants neither TSH nor thyroid disease at 0.05/11 level (10-Feb-2023).
- Rstudio crashed, so that the results are not reproduable anymore (11-Feb-2023).
- We decided to use the natural log-transformed winsorized eGFR (so called eGFRw.log) for interaction analysis rather than __log(eGFR)__. Previously we had found a significant interaction with TSH in one single locus.
```R
# Interaction of SNP-TSH
lm(eGFRw.log ~ SNP * TSH + PC1 + ... + PC10)
# Interaction of SNP-thyroid_disease
lm(eGFRw.log ~ SNP * TSH_cat + PC1 + ... + PC10)
```
- After changing the outcome, we realized if we use as __eGFRw.log__ the outcome, we will have two significant n=interaction terms for the SNPs with TSH in one locus (13-Feb-2023).

- Here we tried to depict the interaction of the dosage with TSH levels for the two variants with significant interaction effect (15-Feb-23).

- Interaction plots for the two variants have been drawn by the interactions package in R, trying to clearly showing the interaction of eGFR-variants with TSH levels by enabling to adjust inside the linear model for age, sex, nad PCs (20-Feb-23).



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

- Fo look-up of the variants we had several options which covered difference genetic databases. These optios are:
1. [Phenoscanner](http://www.phenoscanner.medschl.cam.ac.uk/)
2. [OpenTargets](https://www.opentargets.org/)
3. [GWAS Catalog](https://www.ebi.ac.uk/gwas/home)
4. [PheWeb](https://github.com/statgen/pheweb)

- Having pooled the results of look-up of the variants in phenoscanner with the scan of the variants in OpenTargets provided by David's aid, we found that there is no common traits that remained after interrogating variants association with different traits between these two look-up results and we are still missing the important info like effect size and standard errors of the interrogated variants significant at 5E-08 level. Therefore, we decided to look the variants up in the GWAS catalog database. 

- To do the interrogation of the variants in GWAS Catalog, there are two options: 
(1) is using the [GWASrapidd v. 0.99.14](https://rmagno.eu/gwasrapidd/reference/get_variants.html) which is an R package available on [Cran](https://cran.r-project.org/web/packages/gwasrapidd/) lokking for the variants to find an evidence in GWAS Catalog database.
(2) is the [GWAS Catalog official](https://www.ebi.ac.uk/gwas/docs/api) [REST API](https://github.com/EBISPOT/goci-rest) in Python.

- Here we are taking the advantage of the GWAS Catalog API to automatically query the SNPs to find the corresponding associated traits (14-Feb-2023).

- The summary results of interrogation in GWAS catalog have been obtained and merged with phenoscanner look-up summary results (19-Feb-23).

- Then the joined interrogation summary results was merged with step 3 of the mediation analysis (27-Feb-23).

- The tables in the paper as well as the supplementary files and tables were all updated with the latest finding after merging interrogation results with the mediation analysis steps. The paper was sent to teh PI and pending for his reaction. Then, it would be ready to be circulated among the coauthors (09-Mar-23).

- Thanks to the PI, the paper was circulated among the coauthors (14-Feb-23).

- The feedback from the coauthors were integrated and sent back to the PI. All figures were modified and corrected accordingly (22-Mar-23).

- The last things was to update naming of aPTT in the figure. Then, the paper was submitted to the CHRIS access committee and checked later for the pligerism using LUMC iThenticate tool (28-Mar-23).

# 13. Publication and manuscript revision

- Request for the revison of the paper was received on November 20th.

- List of in-LD variants created using `05-2_linkage_analysis_result.R` (25-Nov-23).

- Extraction of the dosage of in-LD variants was done in five minutes on the Eurac servers using `05-3_extrating_in-LD_variants.sh` (Mon, 19:05, 27-Nov-23).

- Plotted the variance of i-LD variants explained by each PCs `05-4_pca_in-LD_variants.R`. Still need to plot the number of PCs explaining different portion of the total variance of the data (Wed, 18:34, 29-Nov-23).

- Histogram and cumulative distribution of number of PCs explaining different proportions of the total variance of the in-LD variants was genrated and the report was written and sent to the PI (Thu, 20:04, 30-Nov-23).

- The scripts written to be used for running GWAS on ln(eGFRcrea) and doing replication study in the CHRIS study were added to Github (Sat, 23:55, 02-Dec-23).

Dariush

