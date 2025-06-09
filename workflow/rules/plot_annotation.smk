
rule plot_annotation:
	input:
		script = "workflow/scripts/01-6_plot_gennotation.R",
		annotation = ws_path("annotation/{locus}_summary.tsv")
	output:
		plot = ws_path("plot_annotation/{locus}.png")
	conda:
		"../envs/environment.yml"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	params:
		region = "{locus}"
	shell:
		"""
		Rscript {input.script} --annotation {input.annotation} --locus {params.region} --output {output.plot}
		"""
