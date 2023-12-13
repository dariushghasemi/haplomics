

#------------------------#
# list of the desired loci
loci = "SLC34A1 PDILT IGF1R".split()

#------------------------#
# a pseudo-rule that collects the target files
rule all:
	input:
		expand("genotype/{locus}.vcf.gz", locus = loci),
		expand("genotype/{locus}_dosage.txt", locus = loci),
		expand("annotation/{locus}_variants.list", locus = loci),
		expand("annotation/{locus}_annotation.txt", locus = loci),
		expand("output/26-Oct-23_plot_histo_{locus}.png", locus = loci)

#------------------------#
rule get_locus:
	input:
		script = "01-1_subsetting_vcf.sh",
		position = "locus_window.txt"
	output:
		vcf = "genotype/{locus}.vcf.gz",
		sentinel = "annotation/{locus}.sentinel"
	params:
		file = "locus_window.txt"
	shell:
		"""
		bash {input.script} {input.position}
		touch {output.sentinel}
		"""
#------------------------#
rule get_dosage:
	input:
		script = "01-2_extracting_dosage.sh",
		vcf = "genotype/{locus}.vcf.gz",
		sentinel = "annotation/{locus}.sentinel"
	output:
		dosage = "genotype/{locus}_dosage.txt",
		variants = "annotation/{locus}_variants.list",
		annotation = "annotation/{locus}_annotation.txt"
	params: 
		file = "genotype/{locus}.vcf.gz"
	shell:
		"""
		bash {input.script} {input.vcf}
		"""

#------------------------#
rule plot_histogram:
	input:
		script = "01-3_plot_histogram.R",
		variants = "annotation/{locus}_variants.list",
	output:
		plot = "output/26-Oct-23_plot_histo_{locus}.png",

	params:
		file = "annotation/{locus}_variants.list",
	shell:
		"""
		Rscript {input.script} {input.variants}

		"""
#------------------------#
rule build_haplotypes:
	input:
		script = "03-1_building_haplotypes.R",
		dosage = "genotype/{locus}_dosage.txt",
	output:
		result = "output/*_{locus}_*.RDS"
	params:
		dosage = "genotype/{locus}_dosage.txt",
	shell:
		"""
		Rscript {input.script} {input.dosage}
		"""
#------------------------#

ruleorder: get_locus > get_dosage > plot_histogram > build_haplotypes

#------------------------#

