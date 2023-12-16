

#------------------------#
# list of the desired loci
loci = "SLC34A1 PDILT IGF1R".split()

#------------------------#
from datetime import date

today = date.today()
formatted_date = today.strftime('%d-%b-%y')
#print("Today is:", formatted_date)
myday = '15-Dec-23'
#------------------------#
# a pseudo-rule that collects the target files
rule all:
	input:
		expand("genotype/{locus}.vcf.gz", locus = loci),
		expand("genotype/{locus}_dosage.txt", locus = loci),
		expand("annotation/{locus}_variants.list", locus = loci),
		expand("annotation/{locus}_annotation.txt", locus = loci),
		expand("output/26-Oct-23_plot_histo_{locus}.png", locus = loci),
		expand("output/plot_annotations/{day}_{locus}_plot_annotations.png", locus = loci, day = formatted_date),
		expand("output/result_associations/15-Dec-23_{locus}_association_results_full1.RDS", locus = loci),
		expand("output/result_associations/15-Dec-23_{locus}_association_results_short1.RDS", locus = loci),
		expand("output/plot_haplotypes/{day}_{locus}_plot_haplotypes.png", locus = loci, day = formatted_date),
		expand("output/plot_haplotypes/{day}_{locus}_plot_haplotypes_shrinked.png", locus = loci, day = formatted_date),
		expand("output/plot_heatmaps/{day}_{locus}_plot_heatmap_haplotypes_effect.png", locus = loci, day = formatted_date)


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
rule plot_annotation:
	input:
		script = "01-4_plot_annotation.R",
		annotation = "annotation/{locus}_annotation.txt",
	output:
		plot = "output/plot_annotations/{day}_{locus}_plot_annotations.png",
	params:
		file = "annotation/{locus}_annotation.txt",
	shell:
		"""
		Rscript {input.script} {input.annotation}
		"""
#------------------------#
rule build_haplotypes:
	input:
		script = "03-1_haplotypes_building.R",
		dosage = "genotype/{locus}_dosage.txt",
	output:
		result_full  = "output/result_associations/15-Dec-23_{locus}_association_results_full1.RDS",
		result_short = "output/result_associations/15-Dec-23_{locus}_association_results_short1.RDS",
	params:
		file = "genotype/{locus}_dosage.txt",
	shell:
		"""
		Rscript {input.script} {input.dosage}
		"""
#------------------------#
rule plot_haplotypes:
	input:
		script = "03-2_haplotypes_plot.R",
		associations = "output/result_associations/15-Dec-23_{locus}_association_results_short1.RDS",
		annotation = "annotation/{locus}_annotation.txt",
		variants = "annotation/{locus}_variants.list",
	output:
		plot1 = "output/plot_haplotypes/{day}_{locus}_plot_haplotypes.png",
		plot2 = "output/plot_haplotypes/{day}_{locus}_plot_haplotypes_shrinked.png",
	params:
		associations = "output/result_associations/15-Dec-23_{locus}_association_results_short1.RDS",
		annotation = "annotation/{locus}_annotation.txt",
		variants = "annotation/{locus}_variants.list",
	shell:
		"""
		Rscript {input.script} {input.associations} {input.annotation} {input.variants}
		"""
#------------------------#
rule plot_associations:
	input:
		script = "03-3_haplotypes_heatmap.R",
		associations = "output/result_associations/15-Dec-23_{locus}_association_results_full1.RDS",
	output:
		plot = "output/plot_heatmaps/{day}_{locus}_plot_heatmap_haplotypes_effect.png",
	params:
		associations = "output/result_associations/15-Dec-23_{locus}_association_results_full1.RDS",
	shell:
		"""
		Rscript {input.script} {input.associations}
		"""
#------------------------#
ruleorder: get_locus > get_dosage > plot_histogram > plot_annotation > build_haplotypes > plot_haplotypes > plot_associations

#------------------------#

