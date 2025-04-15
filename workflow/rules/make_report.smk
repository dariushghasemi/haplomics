
import os
from os.path import join as pjoin

configfile: "../config/configuration.yml"

# inputs
loci   = config["gene_set"].split()
assays = config["datasets"].split()
resdir = config["path_rep"]
#loci  = "GAB2"


rule all:
	input:
		expand(pjoin(resdir, "results/report_html/{locus}_report.nb.html"), locus = loci)

#------------------------#
rule render_report:
	input:
		markdown = "scripts/04-0_report.Rmd",
		plt_zoom = lambda wildcards: ["/home/dghasemisemeskandeh/projects/gwas/05_regional_association/LZ_plots/09-Mar-23_{locus}-1.png".format(locus=wildcards.locus)],
		plt_hist = lambda wildcards: [pjoin(resdir, "results/plot_histogram/{locus}_plot_histo.png".format(locus=wildcards.locus))],
		plt_anot = lambda wildcards: [pjoin(resdir, "results/plot_annotation/{locus}_plot_annotations.png".format(locus=wildcards.locus))],
		plt_gene = lambda wildcards: [pjoin(resdir, "results/plot_gene/{locus}_plot_gene.png".format(locus=wildcards.locus))],
		plt_hap1 = lambda wildcards: [pjoin(resdir, "results/plot_haplotypes/{locus}/{locus}_{assay}_plot_haplotypes.png".format(locus=wildcards.locus, assay=assay)) for assay in assays],
		plt_hap2 = lambda wildcards: [pjoin(resdir, "results/plot_haplotypes/{locus}/{locus}_{assay}_plot_haplotypes_shrinked.png".format(locus=wildcards.locus, assay=assay)) for assay in assays],
		plt_heat = lambda wildcards: [pjoin(resdir, "results/plot_heatmaps/{locus}/{locus}_{assay}_plot_heatmap.png".format(locus=wildcards.locus, assay=assay)) for assay in assays],
		res_rds  = lambda wildcards: [pjoin(resdir, "results/result_tidied/{locus}/{locus}_{assay}_association_results_tidied.RDS".format(locus=wildcards.locus, assay=assay)) for assay in assays],
		tbl_summ = lambda wildcards: [pjoin(resdir, "results/report/{locus}_merged_data_summary.txt".format(locus=wildcards.locus))],
	output:
		html = pjoin(resdir, "results/report_html/{locus}_report.nb.html")
	params:
		locus = "{locus}",
		assay = assays,
		html = pjoin(resdir, "results/report_html/{locus}_report.nb.html")
	threads: 1
	resources:
		runtime=lambda wc, attempt: attempt * 30,
		mem_mb=4048, disk_mb=20000
	shell:
		"""
		Rscript -e  \
		  'rmarkdown::render(  \
		    input = "{input.markdown}",  \
		    output_file = "{params.html}",  \
		    params = list(  \
		      LOCUS = "{params.locus}",   \
			  ASSAY = "{params.assay}",   \
			  zoom = "{input.plt_zoom}",  \
			  hist = "{input.plt_hist}",  \
			  anot = "{input.plt_anot}",  \
			  gene = "{input.plt_gene}",  \
			  hap1 = "{input.plt_hap1}",  \
			  hap2 = "{input.plt_hap2}",  \
			  heat = "{input.plt_heat}",  \
			  res  = "{input.res_rds}",   \
			  summ = "{input.tbl_summ}" 
		    ))'
		"""
#------------------------#