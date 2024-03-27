#!/usr/bin/bash

#------------------------#
# genetic data originally was built in GRCh37
# but newly version by Michele F is in GRCh38
# so no need for lifting!

#------------------------#
# directories
imputed=/scratch/compgen/data/genetics/CHRIS13K/Imputation/TOPMed/
main=/home/dghasemisemeskandeh/projects/haploAnalysis
window=$1
VCF=$2
odir=$3

# $1 is locus_window.txt
#------------------------#
# Before extracting the variants, 2 steps have been done on the lifted vcf file (1. sorting, 2. indexing)
# 1. bcftools sort $WES/CHRIS13K.WES.imputed.rsq03.vep.hg38.zip.vcf.gz -Oz -o CHRIS13K.WES.imputed.rsq03.vep.hg38.sorted.bgzip.vcf.gz
# 2. tabix CHRIS13K.WES.imputed.rsq03.vep.hg38.sorted.bgzip.vcf.gz

# Then, it was supposed to remove conting from header lines due to lifting VCF, yet dropping conting has consequences! So, I decided to skip it.
#bcftools annotate --remove '##contig=<ID=/d' -Oz9 -o $output_dir/${locus}.vcf.gz $output_dir/${locus}_conting.vcf.gz

#mkdir $output_dir

# remove old files
#rm -r $output_dir/*

#------------------------#
# reading the file with info about window size!
tail -n+2 $window | while IFS=$'\t' read -r chr beg end locus;
    do
    # extract the chrom:beg-end coordinates
    window_size="$chr:$beg"-"$end"
    
    # extract the region from the VCF
    bcftools view $VCF -r ${window_size} -Oz -o $odir &&  \
    #sleep 30

    # creating index file
    tabix $odir

    # Create sentinel file with the correct name
    #touch $output_dir/${locus}.sentinel
done
#------------------------#