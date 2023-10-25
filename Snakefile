
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


rule get_locus:
	input:
		script = "01_subsetting_vcf.sh",
		position = "locus_window.txt"
	output:
		file = "genotype/{locus}.vcf.gz"
	params:
		file = "locus_window.txt"
	shell:
		"""
		bash {input.script} {input.position}
		"""

rule get_dosage:
	input:
		script = "03_extracting_dosage.sh",
		file = "genotype/{locus}.vcf.gz"
	output:
		file = "genotype/{locus}_dosage.txt",
		snps = "genotype/{locus}_variants.list"
	params: 
		file = "genotype/{locus}.vcf.gz"
	shell:
		"""
		bash {input.script} {input.file}
		"""
