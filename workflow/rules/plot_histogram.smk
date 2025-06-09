
rule plot_histogram:
	input:
		script   = "workflow/scripts/01-3_plot_histogram.R",
		variants = ws_path("dosage/{locus}.bim")
	output:
		plot = ws_path("plot_histogram/{locus}.png")
	conda:
		"../envs/environment.yml"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		Rscript {input.script}  \
			--variants {input.variants}  \
			--output {output.plot}
		"""
