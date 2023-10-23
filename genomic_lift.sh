#!/usr/bin/bash

#SBATCH --job-name=lift

chain=/home/demmert/work/Functionalization/work
input=/shared/statgen/CHRIS13K/Imputation/WES
output=/scratch/compgen/data/genetics/CHRIS13K/Imputation/WES

echo CrossMap.py vcf $chain/hg19ToHg38.over.chain.gz \
    $input/CHRIS13K.WES.imputed.rsq03.vep.vcf.gz  \
	--compress	$output/CHRIS13K.WES.imputed.rsq03.vep.hg38.vcf.gz
