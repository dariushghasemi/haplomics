
rule plot_heatmap:
	input:
		script = "workflow/scripts/03-4_haplotypes_heatmap.R",
		result = ws_path("result_associations/{locus}_{dataset}_association_results.RDS")
	output:
		plt = ws_path("plot_heatmaps/{locus}_{dataset}_plot_heatmap.png"),
		tbl = ws_path("result_signif/{locus}_{dataset}_association_results_signif.csv"),
		rds = ws_path("result_tidied/{locus}_{dataset}_association_results_tidied.RDS")
	conda:
		"../envs/environment.yml"
	log:
		ws_path("logs/plot_heatmaps/{locus}_{dataset}.log"),
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		Rscript {input.script} {input.result} {output.plt} {output.tbl} {output.rds} 2> {log}
		"""
