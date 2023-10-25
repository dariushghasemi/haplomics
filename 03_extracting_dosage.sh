#!/usr/bin/bash

# take the vcf file as input
locus=$1

# extract WES variants in the region after merging
bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%DS\n]' ${locus} -o ${locus}_dosage.txt

# extract WES variants dosage after merging with genomAD and gwas 
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' ${locus} -o ${locus}_variants.list

