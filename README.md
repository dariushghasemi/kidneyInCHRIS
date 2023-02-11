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

- We replaced log(eGFR) with eGFRw.log.Res in the above model and ran it again. The results stay consistent in a way that those two variants with significant inetraction with TSH remain significant. At the end, in this way we can avoid over-adjustment of the covariates: age and sex (03-Feb-23).

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

