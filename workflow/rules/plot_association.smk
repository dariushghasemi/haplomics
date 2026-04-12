
rule plot_heatmap:
	input:
		RDS = rules.build_haplotypes.output.result
	output:
		plt = ws_path("plot_heatmaps/{locus}_{dataset}_plot_heatmap.png"),
		rds = ws_path("results/{locus}_{dataset}_tidied.RDS")
	conda:
		"../envs/environment.yml"
	log:
		ws_path("logs/plot_heatmaps/{locus}_{dataset}.log"),
	resources:
		runtime=lambda wc, attempt: attempt * 270,
	script: 
		"../scripts/03-4_haplotypes_heatmap.R"
