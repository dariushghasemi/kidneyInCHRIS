#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import pandas as pd
import sys
import json
import os
import logging
from concurrent.futures import ThreadPoolExecutor
import argparse
from tqdm import tqdm
import time

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Define the GTEx REST API endpoint (updated for API v2)
API_URL = "https://gtexportal.org/api/v2"  # Updated from /rest/v2 to /api/v2
MAX_RETRIES = 3
RETRY_DELAY = 2  # seconds
SIGNIFICANCE_THRESHOLD = 0.05  # Default p-value threshold for significance
MAX_WORKERS = 5  # Limit concurrent requests to avoid overwhelming the API

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Query GTEx API for eQTLs')
    parser.add_argument('--variants_file', help='File containing variant IDs (one per line)')
    parser.add_argument('--output_prefix', help='Prefix for output files')
    parser.add_argument('--dataset', default='gtex_v8', help='GTEx dataset ID (default: gtex_v8)')
    parser.add_argument('--pvalue', type=float, default=SIGNIFICANCE_THRESHOLD, 
                      help=f'P-value threshold for significance (default: {SIGNIFICANCE_THRESHOLD})')
    parser.add_argument('--threads', type=int, default=MAX_WORKERS,
                      help=f'Number of concurrent API requests (default: {MAX_WORKERS})')
    parser.add_argument('--rs-lookup', action='store_true', 
                      help='Input file contains RS IDs that need to be converted to variant IDs')
    return parser.parse_args()

def lookup_variant_id(rs_id):
    """
    Look up variant ID from RS ID using the GTEx variant endpoint.
    Returns the first variant ID found or None if not found.
    """
    url = f"{API_URL}/variant/rsId/{rs_id}"
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200 and response.text:
            data = response.json()
            if data and 'variantId' in data:
                return data['variantId']
    except Exception as e:
        logger.warning(f"Failed to look up variant ID for {rs_id}: {str(e)}")
    return None

def query_variant_eqtl(variant_id, dataset_id, retry_count=0):
    """Query eQTLs for a single variant and return the data."""
    url = f"{API_URL}/association/singleTissueEqtl"
    params = {
        "variantId": variant_id,
        "datasetId": dataset_id
    }
    
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()  # Raise exception for HTTP errors
        
        # Check if the response is valid JSON
        if response.text.strip():
            return variant_id, response.json()
        else:
            raise ValueError("Empty response received")
            
    except requests.exceptions.RequestException as e:
        if retry_count < MAX_RETRIES:
            logger.warning(f"Retrying {variant_id} after error: {str(e)} (attempt {retry_count + 1}/{MAX_RETRIES})")
            time.sleep(RETRY_DELAY * (retry_count + 1))  # Exponential backoff
            return query_variant_eqtl(variant_id, dataset_id, retry_count + 1)
        else:
            logger.error(f"Failed to fetch data for variant {variant_id} after {MAX_RETRIES} attempts: {str(e)}")
            return variant_id, None
    except ValueError as e:
        if retry_count < MAX_RETRIES:
            logger.warning(f"Retrying {variant_id} after error: {str(e)} (attempt {retry_count + 1}/{MAX_RETRIES})")
            time.sleep(RETRY_DELAY * (retry_count + 1))
            return query_variant_eqtl(variant_id, dataset_id, retry_count + 1)
        else:
            logger.error(f"Failed to process response for variant {variant_id} after {MAX_RETRIES} attempts: {str(e)}")
            return variant_id, None

def process_variant_data(variant_data, variant_id, p_value_threshold, original_id=None):
    """Process the data for a variant and return a list of result dictionaries."""
    results = []
    
    if not variant_data or 'data' not in variant_data:
        return results
    
    variant_info = {}
    if original_id:
        variant_info["Original_ID"] = original_id
    
    for result in variant_data['data']:
        # Only include significant associations based on p-value threshold
        p_value = result.get('pValue', 1.0)
        if p_value <= p_value_threshold:
            result_dict = {
                "Variant_id": variant_id,
                "rs_id": result.get('snpId', ''),
                "Dataset": result.get('datasetId', ''),
                "CHR": result.get('chromosome', ''),
                "POS": result.get('pos', ''),
                "Gene_id": result.get('gencodeId', ''),
                "Gene": result.get('geneSymbol', ''),
                "NES": result.get('nes', ''),  # Normalized effect size
                "P_Value": p_value,
                "Tissue": result.get('tissueSiteDetailId', '')
            }
            
            # Add original ID if available
            if original_id:
                result_dict["Original_ID"] = original_id
                
            results.append(result_dict)
    
    return results

def main():
    """Main function to process variants and generate output files."""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Read variants from the text file
    try:
        with open(args.variants_file, "r") as file:
            variants = [line.strip() for line in file if line.strip()]
        
        logger.info(f"Loaded {len(variants)} variants from {args.variants_file}")
    except Exception as e:
        logger.error(f"Error reading variants file: {str(e)}")
        sys.exit(1)
    
    if not variants:
        logger.error("No variants found in the input file")
        sys.exit(1)
    
    # If input is RS IDs, convert them to variant IDs
    variant_mapping = {}
    if args.rs_lookup:
        logger.info("Converting RS IDs to variant IDs...")
        for rs_id in tqdm(variants, desc="Looking up variant IDs"):
            variant_id = lookup_variant_id(rs_id)
            if variant_id:
                variant_mapping[rs_id] = variant_id
                logger.info(f"Mapped {rs_id} to {variant_id}")
            else:
                logger.warning(f"Could not find variant ID for {rs_id}")
        
        # Use the found variant IDs
        if variant_mapping:
            logger.info(f"Successfully mapped {len(variant_mapping)} out of {len(variants)} RS IDs")
            variants_to_query = list(variant_mapping.values())
        else:
            logger.error("Could not map any RS IDs to variant IDs")
            sys.exit(1)
    else:
        # Use the variant IDs directly
        variants_to_query = variants
        
    # Initialize JSON output
    json_output = {}
    
    # Query variants using thread pool for parallel processing
    all_results = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(query_variant_eqtl, variant_id, args.dataset): variant_id for variant_id in variants_to_query}
        
        for future in tqdm(futures, desc="Querying variants"):
            variant_id = futures[future]
            variant_id, data = future.result()
            
            # Find original ID if it was mapped
            original_id = None
            if args.rs_lookup:
                for rs_id, mapped_id in variant_mapping.items():
                    if mapped_id == variant_id:
                        original_id = rs_id
                        break
            
            if data:
                json_output[variant_id] = data
                variant_results = process_variant_data(data, variant_id, args.pvalue, original_id)
                all_results.extend(variant_results)
                if variant_results:
                    logger.info(f"Found {len(variant_results)} significant eQTLs for {variant_id}")
                else:
                    logger.info(f"No significant eQTLs found for {variant_id}")
    
    # Save JSON output
    json_output_file = f"{args.output_prefix}_raw.json"
    with open(json_output_file, "w") as json_file:
        json.dump(json_output, json_file, indent=2)
    
    # Create DataFrame from results
    result_df = pd.DataFrame(all_results)
    
    # Sort results by p-value (most significant first)
    if not result_df.empty:
        result_df = result_df.sort_values(by=['P_Value', 'Gene', 'Tissue'])
    
    # Save summary to TSV
    tsv_output_file = f"{args.output_prefix}_summary.tsv"
    result_df.to_csv(tsv_output_file, sep="\t", index=False)
    
    # Save tissue-specific summaries
    if not result_df.empty:
        # Count significant associations per tissue
        tissue_counts = result_df['Tissue'].value_counts().reset_index()
        tissue_counts.columns = ['Tissue', 'Significant_eQTLs']
        tissue_output_file = f"{args.output_prefix}_tissue_summary.tsv"
        tissue_counts.to_csv(tissue_output_file, sep="\t", index=False)
        
        # Count significant associations per gene
        gene_counts = result_df['Gene'].value_counts().reset_index()
        gene_counts.columns = ['Gene', 'Significant_eQTLs']
        gene_output_file = f"{args.output_prefix}_gene_summary.tsv"
        gene_counts.to_csv(gene_output_file, sep="\t", index=False)
    
    # Print summary statistics
    logger.info(f"\nGTEx eQTL query complete!")
    logger.info(f"Found {len(all_results)} significant eQTL associations (p ≤ {args.pvalue})")
    if not result_df.empty:
        logger.info(f"Unique genes: {result_df['Gene'].nunique()}")
        logger.info(f"Unique tissues: {result_df['Tissue'].nunique()}")
    
    logger.info(f"\nOutput files:")
    logger.info(f"  Raw JSON data: {json_output_file}")
    logger.info(f"  Summary TSV: {tsv_output_file}")
    if not result_df.empty:
        logger.info(f"  Tissue summary: {tissue_output_file}")
        logger.info(f"  Gene summary: {gene_output_file}")

if __name__ == "__main__":
    main()