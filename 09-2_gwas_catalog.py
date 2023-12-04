#!/usr/bin/env python3

import requests
import json
import pandas as pd
import numpy as np
import time
import csv
from pandas import json_normalize

#print(response.status_code)
#print(response.headers)
#print(response.headers["Content-Type"])
#----------------#

# Define function to retrieve info from teh rest API
def query_gwas(variants):
    summary = []

    for variant in variants:
        url = "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/" + variant + "/associations/"
        # JSON output of the url
        response = requests.get(url)

        if response.status_code == 200:
            # Access the data in the JSON
            # and parse the JSON output
            associations        = json.loads(response.text)['_embedded']['associations']
            for association in associations:
            
                if len(association['loci'][0]['authorReportedGenes']) > 0:
                
                    locus       = association['loci'][0]['authorReportedGenes'][0]['geneName']
                    risk_allele = association['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName']
                    risk_A      = risk_allele.split('-')[1]
                else:
                    "N/A"
                # Access fields as needed
                risk_AF         = association['riskFrequency']
                beta            = association['betaNum']
                stderr          = association['standardError']
                pvalue          = association['pvalue']
                study_url       = association['_links']['study']['href']
                # Retrieve the study trait from the study URL
                study_response  = requests.get(study_url)
                if study_response.status_code == 200:
                
                    study_json  = json.loads(study_response.text)
                    study_trait = study_json['diseaseTrait']['trait']
                    ancestry    = study_json['ancestries'][0]['ancestralGroups'][0]['ancestralGroup']
                    pubmedId    = study_json['publicationInfo']['pubmedId']
                    
                else:
                    study_trait = ""
                #print("Traits:", ", ".join([d["trait"] for d in data]))
                summary.append([variant, locus, risk_A, risk_AF, beta, stderr, pvalue, study_trait, ancestry, pubmedId])
    # Create the summary as dataframe
    df = pd.DataFrame(summary, columns = ["variant", "Locus", "Risk_A", "Risk_AF", "Beta", "Stderr", "Pvalue", "Trait", "Ancestry", "Pubmed_ID"])
    return df
#----------------#

# Importing the replicated variants rsIDs
repSNPs = pd.read_table("/home/dghasemisemeskandeh/projects/gwas/replicationAnalysis/scripts/output/10-Jan-23_replicatedSNPs_eGFRw.log.Res.txt", delimiter = "\t")

# Define the list of variants of interest
#variants = ["rs7329174", "rs3812036", "rs13146355"]]

# Use df=df.replace(np.nan,0,regex=True) function
# to replace the ‘NaN’ values with ‘0’ value

repSNPs1 = repSNPs.replace(np.nan, "NA", regex = True)
repSNPs2 = repSNPs.dropna(subset = ['RSID'], inplace = False)

mySNPs   = repSNPs2['RSID'].values.tolist() #astype(float)

print("\n| N. o SNPs | ", "N. of SNPs with non-missing rsID |") 
print("|-----------------------------------------------|")
print("|     ", len(repSNPs['RSID']), " | ", len(mySNPs), " |")

#----------------#

# Iterate the function
df = query_gwas(mySNPs)  

# Dimention of the summary results
print(df.size)

# Save the summary results of variants interrogation in GWAS Catalog
df.to_csv(r'~/projects/variant_interrogation/24-Feb-23_GWAS Catalog interrogation summary results.csv', index = False, header = True) # , sep = "\t"

print("\n\n Note: \n", len(np.unique(df['variant'])), "variants", "out of", len(mySNPs), "were successfully found in the databse")
#----------------#
