# haploAnalysis
Haplotype analysis pipeline for the replicated genetic kidney loci in CHRIS study

- Attaining the total number of variants (Sat, 17:42, 22-Oct-23):
 
```bash
# imputed WES - b37
bcftools index -n /scratch/compgen/data/genetics/CHRIS13K/Imputation/WES/CHRIS13K.WES.imputed.rsq03.vep.vcf.gz
1034420

# lifted imputed WES - b38
bcftools index -n CHRIS13K.WES.imputed.rsq03.vep.hg38.sorted.bgzip.vcf.gz
1033325

# lifted and sorted imputed WES - b38
zcat /scratch/compgen/data/genetics/CHRIS13K/Imputation/WES/CHRIS13K.WES.imputed.rsq03.vep.hg38.zip.vcf.gz | cut -f1-3 | grep -v "#" | wc -l
1033325
```

- activate conda and snakemake environmemnt

```bash
export MAMBA_ROOT_PREFIX=/shared/Software/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate snakemake

# to see the activated environments
micromamba env list
```
# Analysis steps
- Step 1: defining the window for each locus respect to the recombination rate
- Step 2: extracting the variants within the window (or recombination spikes)
- Step 3: taking the dosage level of the extracted variants for each locus
- Step 4: 
- Step 5: 

```bash
rule targets:
	input:
		get_ranges()

def get_ranges(wildcards):
	myfiles = "genotype/{locus}.vcf.gz"

	return(myfiles)

rule two:
	input:
		get_ranges()


direct = directory("genotype/{locus}"),

output:
	odir = directory("genotype/{locus}"),
	ofile = "genotype/{locus}/{locus}.vcf.gz"

rule all:
	input:
		"genotype/PDILT.vcf.gz",
		"genotype/SLC34A1.vcf.gz",
		"genotype/IGF1R.vcf.gz",
		"genotype/PDILT_dosage.txt",
		"genotype/SLC34A1_dosage.txt",
		"genotype/IGF1R_dosage.txt",
		"genotype/PDILT_variants.list",
		"genotype/SLC34A1_variants.list",
		"genotype/IGF1R_variants.list"

# or 
# Determine the {vcf} for all loci
#vcf, = glob_wildcards("genotype/{locus}.vcf.gz")
#dosag, = glob_wildcards("genotype/{locus}_dosage.txt")
#annot, = glob_wildcards("genotype/{locus}_variants.list")

```
- Allele frequency of the extracted variants at each locus was depicted by a histogram (Thu, 18:42, 26-Oct-23).