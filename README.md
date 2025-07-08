<img src="https://github.com/DariushG3/kidneyInCHRIS/blob/main/Untitled.jpg" width=250 align="right">

## KidneyInCHRIS: A Part of TrainCKDis Horizon 2020 EU Project

### Genetic Determinants of Kidney Function in the CHRIS Study
This project, started in February 2021, investigates the genetic determinants of kidney function using data from the **C**orporative **H**ealth **R**esearch **I**n **S**outh Tyrol (CHRIS) study. Our research focuses on genome-wide association studies (GWAS) and post-GWAS statistical analyses to explore genetic variations influencing kidney function.
___________________________________________________________________________________________________________

#### Project Overview
We analyze 147 loci associated with creatinine-based eGFR, identified by CKDGen’s European ancestry GWAS meta-analysis ([Wuttke et al., Nat Genet 2019](https://www.nature.com/articles/s41588-019-0407-x)). Our key objectives include:

- Identifying significant genetic loci relevant to kidney function.
- Conducting GWAS on renal function traits.
- Comparing effect sizes of replicated variants with discovery studies.
- Investigating interactions between kidney function and thyroid-related traits.
- Performing sensitivity analyses to assess robustness of findings.

#### Methodology
The methods are thoroughly described in the paper. To summarize the methods, we followed several analyses steps as follows:
1. **Selection of Genetic Loci**: We selected lead variants from 147 loci identified in CKDGen and retained proxy SNPs (LD r² > 0.8) within a 1M bp region around each lead variant in CHRIS.
2. **GWAS on Renal Function**: We performed a genome-wide association study (GWAS) mapping CHRIS 10K participants’ genotypes to the phenotype ln(eGFR).
3. **Comparison of Replicated Variants**: After comparing the loci in CHRIS 10K participants and in the CKDGen using regional association plot, we evaluated the magnitude of effect sizes for replicated variants compared to discovery studies, assessing their consistency and significance respect of minor allele frequency differences.
4. **Phenome-Wide Mediation Analysis**: We implemented a comprehensive phenome-wide mediation analysis to assess how other biological pathways influence kidney function.
5. **Kidney-Thyroid Interaction Analysis**: We explored interactions between kidney function and thyroid-related traits (TSH/FT3/FT4 levels) to investigate potential biological mechanisms.
6. **Sensitivity Analyses**: To ensure robustness, we conducted sensitivity tests:
    - SNP Dosage Sensitivity to FT3/FT4 Missingness: We assessed whether missing FT3/FT4 data influenced SNP dosage levels using a Wilcoxon's test.
    - Municipality-Adjusted SNP-Kidney Associations: To account for regional differences, we adjusted SNP-kidney associations for municipality to detect potential batch effects.
7. **Interrogation of SNPs in Phenoscanner**: We cross-referenced 163 replicated SNPs with Phenoscanner to explore their associations with other traits and diseases.

___________________________________________________________________________________________________________
#### Publication & Future Work
The findings from this project contribute to the ongoing TrainCKDis Horizon 2020 EU Project. Our work is currently published on medrxiv with further refinements based on peer review and manuscript revisions:
- **[preprinted manuscript version 1](https://www.medrxiv.org/content/10.1101/2023.04.15.23288540v1)**
- **[preprinted manuscript version 2](https://www.medrxiv.org/content/10.1101/2023.04.15.23288540v2)**
- **[published manuscript on PLoS One](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0323057)**


Dariush Ghasemi
