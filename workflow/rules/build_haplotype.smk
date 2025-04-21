
rule build_haplotypes:
	input:
		script = "workflow/scripts/03-2_haplotypes_building.R",
		data   = "results/merged_data/{locus}_{dataset}_merged_data.RDS"
	output:
		result = "results/result_associations/{locus}_{dataset}_association_results.RDS"
	conda:
		"../envs/environment.yml"
	params:
		covar_file = config["covariates_file"],
		min_freq = config.get("thresholds").get("min_freq"),
		max_haps = config.get("thresholds").get("max_haps"),
		min_pp = config.get("thresholds").get("min_pp"),
		n_batch = 2,
		n_try = 2,
	log:
		"logs/build_haplotypes/{locus}_{dataset}.log"
	resources:
		runtime=lambda wc, attempt: attempt * 6000,
		mem_mb=get_mem_mb
	shell:
		"""
		Rscript {input.script}  \
			--data {input.data}  \
			--covariate {params.covar_file} \
			--min_freq {params.min_freq} \
			--max_haps {params.max_haps} \
			--min_pp {params.min_pp} \
			--n_batch {params.n_batch} \
			--n_try {params.n_try} \
			--output {output.result} 2> {log}
		"""
