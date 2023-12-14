#!/usr/bin/bash

# take the vcf file as input
locus_vcf=$1

# base name of the locus file
locus_file="$(basename -- $locus_vcf)"

#------------------------#
# removing file format
locus=$(echo $locus_file | sed -e 's/genotype.//g' -e 's/.vcf.gz//g')

# attain the number of exonic variants in the locus
n_SNPs=$(bcftools index -n $locus_vcf)

# briefly describing the process in this script
echo Working on "${n_SNPs}" exonic variants at "${locus}" locus...

#------------------------#
# extract WES variants in the region after merging
bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%DS\n]' genotype/${locus}.vcf.gz -o genotype/${locus}_dosage.txt

# extract WES variants dosage after merging with genomAD and gwas 
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' genotype/${locus}.vcf.gz -o annotation/${locus}_variants.list

# extract variants annotation
bcftools +split-vep  -s worst -f '%CHROM\t%POS\t%ID\t%SYMBOL\t%Gene\t%Consequence\n' genotype/${locus}.vcf.gz -o annotation/${locus}_annotation.txt

#------------------------#
# ending message!
echo The dosage levels, characteristics, and annotations of the variants in "${locus}" locus were created.
