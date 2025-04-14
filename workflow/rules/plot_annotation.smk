
rule plot_annotation:
	input:
		script = "scripts/01-4_plot_annotation.R",
		annotation = "results/annotation/{locus}_annotation.txt"
	output:
		plot = "results/plot_annotation/{locus}_plot_annotations.png"
	conda:
		"../envs/environment.yml"
	params:
		#annotation
	shell:
		"""
		Rscript {input.script} {input.annotation} {output.plot}
		"""
