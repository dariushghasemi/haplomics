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
- Allele frequency of the extracted variants at each locus was depicted using `01-3_plot_histogram.R` (Thu, 18:42, 26-Oct-23).

- Add haplotype reconstruction script to the pipeline using `03-1_haplotypes_building.R` (Wed, 16:50, 13-Dec-23).

- Draw the workflow diagram (Wed, 13-Dec-23).

```bash
snakemake --dag targets | dot -Tpng > Tag.png
```

- Add visualization of the recnstructed haplotype to the workflow (Wed, 21:30, 13-Dec-23).

- Haplotypes plot were generated in two version, full variants and shrinked variants, using `03-2_plot_haplotypes.R` (Thu, 20:16, 14-Dec-23)

- Variants consequences were visulalized for each locus using `01-4_plot_annotation.R` (Fri, 02:02, 15-Dec-23).

- Working to add heatmap plot to visualize the haplotypes associations results using `03-3_haplotypes_heatmap.R` (Fri, 18:10, 15-Dec-23).

- Meanwhile, I'm going attend at Eurac X-Mass party at NOI techpark (Fri, 18:10, 15-Dec-23).

- Heatmap plots of haplotypes associations with the traits were added to the pipeline using `03-3_haplotypes_heatmap.R`.

- The Rmarkdown report was generated using `04-0_report.Rmd`.

- The Rmarkdown report was automated for the loci using `04-1_report_run.sh` (Sat, 22:16, 16-Dec-23).

```bash
# limit the job run to just one rule (-n for a dry-run)
snakemake -R --until plot_associations -n

```
- The frequency of the variants of each gene at the input locus was depicted using `01-5_plot_genes.R` (Mon, 14:20, 18-Dec-23).

- Trying to run the jobs on clusters using below (Mon, 18:40, 18-dec-23):
```bash
snakemake --latency-wait 60 --use-conda --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition}  --mem-per-cpu={cluster.mem} -c {cluster.cores}" --jobs 20
```

Dariush