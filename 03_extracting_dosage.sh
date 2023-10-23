#!/usr/bin/bash

# extract WES variants in the region after merging
bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%DS\n]' genotype/${locus}.vcf.gz -o genotype/${locus}_dosage.txt

# extract WES variants dosage after merging with genomAD and gwas 
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' genotype/${locus}.vcf.gz -o genotype/${locus}_variants.list

