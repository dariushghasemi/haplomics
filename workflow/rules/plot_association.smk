
rule plot_heatmap:
	input:
		script = "workflow/scripts/03-4_haplotypes_heatmap.R",
		result = "results/result_associations/{locus}_{dataset}_association_results.RDS"
	output:
		plt = "results/plot_heatmaps/{locus}_{dataset}_plot_heatmap.png",
		tbl = "results/result_signif/{locus}_{dataset}_association_results_signif.csv",
		rds = "results/result_tidied/{locus}_{dataset}_association_results_tidied.RDS"
	conda:
		"../envs/environment.yml"
	params:
		#result
	log:
		"logs/plot_heatmaps/{locus}_{dataset}.log"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
		mem_mb=get_mem_plt, disk_mb=32000
	shell:
		"""
		Rscript {input.script} {input.result} {output.plt} {output.tbl} {output.rds} 2> {log}
		"""
