
rule plot_gene:
	input:
		script = "scripts/01-5_plot_genes.R",
		annotation = "results/annotation/{locus}_annotation.txt"
	output:
		plot = "results/plot_gene/{locus}_plot_gene.png"
	conda:
		"../envs/environment.yml"
	params:
		#annotation
	shell:
		"""
		Rscript {input.script} {input.annotation} {output.plot}
		"""
