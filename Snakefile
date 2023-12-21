
'''
snakemake   --reason  \
			--until get_locus  \
			--jobs 3  \
			--default-resource mem_gb=8GB  \
			--latency-wait 30  \
			--keep-going  \
			--cluster '
				sbatch  --partition fast  \
						--cores 3  \
						--mem-per-cpu=8800 \
						--output  output/{rule}.{wildcards}.out \
						--error   output/{rule}.{wildcards}.err'
'''
#{resources.mem_gb} \

#------------------------#
# list of the desired loci
loci = "SLC34A1 PDILT IGF1R".split()

#------------------------#
from datetime import date

today = date.today()
my_date = today.strftime('%d-%b-%y')
#print("Today is:", formatted_date)

#------------------------#
# a pseudo-rule that collects the target files
rule all:
	input:
		expand("data/genotype/{locus}.vcf.gz", locus = loci),
		expand("data/dosage/{locus}_dosage.txt", locus = loci),
		expand("data/annotation/{locus}_variants.list", locus = loci),
		expand("data/annotation/{locus}_annotation.txt", locus = loci),
		expand("output/plot_histogram/{day}_{locus}_plot_histo.png", locus = loci, day = my_date),
		expand("output/plot_annotations/{day}_{locus}_plot_annotations.png", locus = loci, day = my_date),
		expand("output/plot_genes/{day}_{locus}_plot_genes.png", locus = loci, day = my_date),
		expand("data/{locus}_haplotypes_data.csv", locus = loci, day = my_date),
		expand("output/result_associations/{locus}_association_results.RDS", day = my_date, locus = loci),
		expand("output/plot_haplotypes/{day}_{locus}_plot_haplotypes.png", locus = loci, day = my_date),
		expand("output/plot_haplotypes/{day}_{locus}_plot_haplotypes_shrinked.png", locus = loci, day = my_date),
		expand("output/plot_heatmaps/{day}_{locus}_plot_heatmap_haplotypes_effect.png", locus = loci, day = my_date),
		expand("output/report_html/{day}_{locus}_report.nb.html", locus = loci, day = my_date)


#------------------------#
rule get_locus:
	input:
		script = "01-1_subsetting_vcf.sh",
		position = "data/locus_window.txt"
	output:
		vcf = "data/genotype/{locus}.vcf.gz",
		sentinel = "data/annotation/{locus}.sentinel"
	params:
		position = "data/locus_window.txt"
	shell:
		"""
		sleep 10
		bash {input.script} {input.position}
		touch {output.sentinel}
		"""
#------------------------#
rule get_dosage:
	input:
		script = "01-2_extracting_dosage.sh",
		vcf = "data/genotype/{locus}.vcf.gz",
		sentinel = "data/annotation/{locus}.sentinel"
	output:
		dosage = "data/dosage/{locus}_dosage.txt",
		variants = "data/annotation/{locus}_variants.list",
		annotation = "data/annotation/{locus}_annotation.txt"
	params: 
		vcf = "data/genotype/{locus}.vcf.gz",
		sentinel = "data/annotation/{locus}.sentinel"
	shell:
		"""
		bash {input.script} {input.vcf}
		"""
#------------------------#
rule plot_histogram:
	input:
		script = "01-3_plot_histogram.R",
		variants = "data/annotation/{locus}_variants.list"
	output:
		plot = "output/plot_histogram/{day}_{locus}_plot_histo.png"
	params:
		variants = "data/annotation/{locus}_variants.list"
	shell:
		"""
		Rscript {input.script} {input.variants}
		"""
#------------------------#
rule plot_annotation:
	input:
		script = "01-4_plot_annotation.R",
		annotation = "data/annotation/{locus}_annotation.txt"
	output:
		plot = "output/plot_annotations/{day}_{locus}_plot_annotations.png"
	params:
		annotation = "data/annotation/{locus}_annotation.txt"
	shell:
		"""
		Rscript {input.script} {input.annotation}
		"""
#------------------------#
rule plot_genes:
	input:
		script = "01-5_plot_genes.R",
		annotation = "data/annotation/{locus}_annotation.txt"
	output:
		plot = "output/plot_genes/{day}_{locus}_plot_genes.png"
	params:
		annotation = "data/annotation/{locus}_annotation.txt"
	shell:
		"""
		Rscript {input.script} {input.annotation}
		"""
#------------------------#
rule build_haplotypes:
	input:
		script = "03-1_haplotypes_building.R",
		dosage = "data/dosage/{locus}_dosage.txt"
	output:
		haplo_data   = "data/{locus}_haplotypes_data.csv",
		result       = "output/result_associations/{locus}_association_results.RDS"
	params:
		dosage = "data/dosage/{locus}_dosage.txt"
	shell:
		"""
		Rscript {input.script} {input.dosage}
		"""
#------------------------#
rule plot_haplotypes:
	input:
		script = "03-2_haplotypes_plot.R",
		result = "output/result_associations/{day}_{locus}_association_results.RDS",
		annotation = "data/annotation/{locus}_annotation.txt",
		variants = "data/annotation/{locus}_variants.list"
	output:
		plot1 = "output/plot_haplotypes/{day}_{locus}_plot_haplotypes.png",
		plot2 = "output/plot_haplotypes/{day}_{locus}_plot_haplotypes_shrinked.png"
	params:
		result = "output/result_associations/{day}_{locus}_association_results.RDS",
		annotation = "data/annotation/{locus}_annotation.txt",
		variants = "data/annotation/{locus}_variants.list"
	shell:
		"""
		Rscript {input.script} {input.result} {input.annotation} {input.variants}
		"""
#------------------------#
rule plot_associations:
	input:
		script = "03-3_haplotypes_heatmap.R",
		result = "output/result_associations/{day}_{locus}_association_results.RDS"
	output:
		plot = "output/plot_heatmaps/{day}_{locus}_plot_heatmap_haplotypes_effect.png"
	params:
		result = "output/result_associations/{day}_{locus}_association_results.RDS"
	shell:
		"""
		Rscript {input.script} {input.result}
		"""
#------------------------#
rule render_report:
    input:
        script = "04-1_report_run.sh",
        reprt  = "04-0_report.Rmd"
    output:
        html = "{day}_{locus}_report.nb.html"
    params:
        #locus = "{locus}"
    shell:
        """
		bash {input.script} {params.locus}
		"""
#------------------------#

ruleorder: get_locus > get_dosage > plot_histogram > plot_annotation > plot_genes > build_haplotypes > plot_haplotypes > plot_associations > render_report

#------------------------#

