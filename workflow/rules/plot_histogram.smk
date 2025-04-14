
rule plot_histogram:
	input:
		script   = "scripts/01-3_plot_histogram.R",
		variants = "results/annotation/{locus}_variants.list"
	output:
		plot = "results/plot_histogram/{locus}_plot_histo.png"
	conda:
		"../envs/environment.yml"
	shell:
		"""
		Rscript {input.script} {input.variants} {output.plot}
		"""
