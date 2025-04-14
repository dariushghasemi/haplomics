#!/bin/python

import requests
import json
import time
from typing import List, Dict, Tuple
import pandas as pd
from pathlib import Path

class VEPAnnotator:
    def __init__(self, server: str = "https://rest.ensembl.org"):
        self.server = server
        self.headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }

    def parse_variant(self, variant_str: str) -> Tuple[str, int, str, str]:
        parts = variant_str.strip().split(':')
        if len(parts) != 4:
            raise ValueError(f"Invalid variant format: {variant_str}")
        return parts[0], int(parts[1]), parts[2], parts[3]

    def format_allele_string(self, ref: str, alt: str) -> str:
        """Format allele string for VEP API"""
        # Handle insertions
        if len(alt) > len(ref) and alt.startswith(ref):
            return alt[len(ref):]
        # Handle deletions
        elif len(ref) > len(alt) and ref.startswith(alt):
            return '-'
        # Handle SNPs and other variants
        return alt

    def annotate_variant(self, chrom: str, position: int, ref_allele: str, alt_allele: str) -> Dict:
        """
        Annotate a single variant using VEP API with automatic allele swapping if needed
        """
        variant_id = f"{chrom}:{position}:{ref_allele}:{alt_allele}"
        print(f"Processing variant: {variant_id}")
        
        # First try with alt_allele
        api_allele = self.format_allele_string(ref_allele, alt_allele)
        endpoint = f"/vep/human/region/{chrom}:{position}-{position}/{api_allele}?"
        
        try:
            response = requests.get(
                self.server + endpoint,
                headers=self.headers
            )
            
            # If we get a reference allele error, try swapping the alleles
            if response.status_code == 400 and "matches reference" in response.text:
                print(f"  Alternative allele matches reference, swapping alleles for {variant_id}")
                api_allele = self.format_allele_string(alt_allele, ref_allele)
                endpoint = f"/vep/human/region/{chrom}:{position}-{position}/{api_allele}?"
                
                response = requests.get(
                    self.server + endpoint,
                    headers=self.headers
                )
            
            if response.ok:
                data = response.json()
                if data:
                    data[0].update({
                        'input_ref_allele': ref_allele,
                        'input_alt_allele': alt_allele
                    })
                    print(f"  Successfully annotated variant: {variant_id}")
                    return data
                else:
                    print(f"  No annotation data returned for variant: {variant_id}")
                    return None
            
            if response.status_code == 429:
                retry_after = int(response.headers.get('Retry-After', 60))
                print(f"  Rate limit hit, waiting {retry_after} seconds...")
                time.sleep(retry_after)
                return self.annotate_variant(chrom, position, ref_allele, alt_allele)
            
            print(f"  API Error for {variant_id}: {response.status_code} - {response.text}")
            return None
                
        except requests.exceptions.RequestException as e:
            print(f"  Error annotating variant {variant_id}: {str(e)}")
            return None
            
        time.sleep(1)
        return None

    def read_variants_from_file(self, file_path: str) -> List[str]:
        """Read variants from a text file"""
        try:
            with open(file_path, 'r') as f:
                # Read lines and remove empty lines and whitespace
                variants = [line.strip() for line in f if line.strip()]
            print(f"Read {len(variants)} variants from {file_path}")
            return variants
        except Exception as e:
            print(f"Error reading variants file: {str(e)}")
            return []

    def process_variant_list(self, variants: List[str]) -> List[Dict]:
        """Process a list of variants and return annotations"""
        annotations = []
        total_variants = len(variants)
        
        print(f"\nStarting annotation of {total_variants} variants...")
        for i, variant in enumerate(variants, 1):
            try:
                chrom, pos, ref, alt = self.parse_variant(variant)
                result = self.annotate_variant(chrom, pos, ref, alt)
                
                if result:
                    annotations.append({
                        'variant': variant,
                        'annotation': result
                    })
                    
            except ValueError as e:
                print(f"  Error processing variant {variant}: {str(e)}")
                continue
            
            print(f"Progress: {i}/{total_variants} variants processed\n")
                
        print(f"Successfully annotated {len(annotations)} out of {total_variants} variants")
        return annotations

    def extract_key_information(self, annotations: List[Dict]) -> pd.DataFrame:
        """Extract key information from annotations into a DataFrame"""
        # [Previous implementation remains the same]
        rows = []
        for entry in annotations:
            variant = entry['variant']
            if not entry.get('annotation'):
                continue
                
            anno = entry['annotation'][0]
            
            # Get consequences from all transcripts
            all_consequences = set()
            coding_consequences = set()
            protein_changes = set()
            
            for trans in anno.get('transcript_consequences', []):
                cons = trans.get('consequence_terms', [])
                all_consequences.update(cons)
                
                if trans.get('biotype') == 'protein_coding':
                    coding_consequences.update(cons)
                    
                    # Collect protein changes
                    if trans.get('protein_start'):
                        aa_change = trans.get('amino_acids', '')
                        hgvsp = trans.get('hgvsp', '')
                        if aa_change or hgvsp:
                            protein_changes.add(f"{trans.get('gene_symbol', '')}:"
                                             f"{aa_change if aa_change else hgvsp}")
            
            # Get the first protein-coding transcript for detailed info
            coding_trans = [t for t in anno.get('transcript_consequences', []) 
                          if t.get('biotype') == 'protein_coding']
            main_trans = coding_trans[0] if coding_trans else anno.get('transcript_consequences', [{}])[0]
            
            # Get clinical significance and RS IDs
            clin_sig = []
            rs_ids = []
            for var in anno.get('colocated_variants', []):
                if var.get('clin_sig'):
                    clin_sig.extend(var['clin_sig'])
                if var.get('id') and var['id'].startswith('rs'):
                    rs_ids.append(var['id'])
            
            # Get frequencies
            freqs = {}
            for var in anno.get('colocated_variants', []):
                if var.get('frequencies'):
                    for allele, freq_data in var['frequencies'].items():
                        if allele == anno.get('input_alt_allele'):
                            freqs = freq_data
                            break
            
            row = {
                'variant': variant,
                'rs_ids': ';'.join(rs_ids),
                'ref_allele': anno.get('input_ref_allele', ''),
                'alt_allele': anno.get('input_alt_allele', ''),
                'most_severe_consequence': anno.get('most_severe_consequence', ''),
                'all_consequences': ';'.join(sorted(all_consequences)),
                'coding_consequences': ';'.join(sorted(coding_consequences)),
                'protein_changes': ';'.join(sorted(protein_changes)),
                'impact': main_trans.get('impact', ''),
                'gene_id': main_trans.get('gene_id', ''),
                'gene_symbol': main_trans.get('gene_symbol', ''),
                'transcript_id': main_trans.get('transcript_id', ''),
                'canonical_transcript': 'YES' if main_trans.get('canonical') == 1 else 'NO',
                'biotype': main_trans.get('biotype', ''),
                'hgnc_id': main_trans.get('hgnc_id', ''),
                'clinical_significance': ';'.join(sorted(set(clin_sig))) if clin_sig else '',
                'sift_prediction': main_trans.get('sift_prediction', ''),
                'polyphen_prediction': main_trans.get('polyphen_prediction', ''),
                'gnomad_af': freqs.get('gnomadg', ''),
                'gnomad_afr_af': freqs.get('gnomadg_afr', ''),
                'gnomad_eur_af': freqs.get('gnomadg_nfe', ''),
                'gnomad_eas_af': freqs.get('gnomadg_eas', '')
            }
            rows.append(row)
            
        return pd.DataFrame(rows)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Annotate variants using VEP API')
    parser.add_argument('--snps-list', help='Path to input variants list (one per line in format chr:pos:ref:alt)')
    parser.add_argument('--output-summary', default='NULL', help='Required name for VEP annotation results')
    parser.add_argument('--output-json', default='NULL', help='Required name for VEP outputs')
    
    args = parser.parse_args()
    
    annotator = VEPAnnotator()
    
    # Read variants from file
    variants = annotator.read_variants_from_file(args.snps_list)
    if not variants:
        print("No variants found in input file. Exiting.")
        return
    
    # Process variants
    annotations = annotator.process_variant_list(variants)
    
    # Save results
    json_output = f"{args.output_json}"
    tsv_output = f"{args.output_summary}"
    
    with open(json_output, 'w') as f:
        json.dump(annotations, f, indent=2)
    
    df = annotator.extract_key_information(annotations)
    df.to_csv(tsv_output, sep='\t', index=False)
    
    print(f"\nAnnotation complete!")
    print(f"Full annotations saved to: {json_output}")
    print(f"Summary saved to: {tsv_output}")
    
if __name__ == "__main__":
    main()