
import os
from os.path import join as pjoin


configfile: "../config/configuration.yml"

# inputs
#loci  = config["gene_set1"].split()
loci  = "STC1"
assays = config["datasets"].split()
resdir = config["path_rep"]


rule all:
	input:
		expand(pjoin(resdir, "results/{locus}_report.nb.html"), locus = loci)


#------------------------#
rule render_report:
	input:
		markdown = "scripts/04-0_report.Rmd",
		plt_zoom = expand("/home/dghasemisemeskandeh/projects/gwas/05_regional_association/LZ_plots/09-Mar-23_{locus}-1.png", locus = loci),
		plt_hist = expand(pjoin(resdir, "results/plot_histogram/{locus}_plot_histo.png"), locus = loci),
		plt_anot = expand(pjoin(resdir, "results/plot_annotation/{locus}_plot_annotations.png"), locus = loci),
		plt_gene = expand(pjoin(resdir, "results/plot_genes/{locus}_plot_genes.png"), locus = loci),
		plt_hap1 = expand(pjoin(resdir, "results/old/plot_haplotypes/{locus}_{assay}_plot_haplotypes.png"), locus = loci, assay = assays),
		plt_hap2 = expand(pjoin(resdir, "results/old/plot_haplotypes/{locus}_{assay}_plot_haplotypes_shrinked.png"), locus = loci, assay = assays),
		plt_heat = expand(pjoin(resdir, "results/old/plot_heatmaps/{locus}_{assay}_plot_heatmap_haplotypes_effect.png"), locus = loci, assay = assays),
		res_sig  = expand(pjoin(resdir, "../output/significant_result/{locus}_{assay}_association_results.csv"), locus = loci, assay = assays)
	output:
		html = pjoin(resdir, "results/{locus}_report.nb.html")
	params:
		locus = "{locus}",
		assays= assays,
		html = pjoin(resdir, "results/{locus}_report.nb.html")
	shell:
		"""
		Rscript -e  \
		  'rmarkdown::render(  \
		    input = "{input.markdown}",  \
		    output_file = "{params.html}",  \
		    params = list(  \
		      LOCUS = "{params.locus}",   \
			  ASSAY = "{params.assays}",   \
			  zoom = "{input.plt_zoom}",  \
			  hist = "{input.plt_hist}",  \
			  anot = "{input.plt_anot}",  \
			  gene = "{input.plt_gene}",  \
			  hap1 = "{input.plt_hap1}",  \
			  hap2 = "{input.plt_hap2}",  \
			  heat = "{input.plt_heat}",  \
			  res  = "{input.res_sig}"
		    ))'
		"""
