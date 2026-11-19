
rule execute_glm:
	input:
		prefix = rules.compute_matrix.output.sentinel,
	output:
		result = ws_path("results/{locus}_{dataset}_phased_results.RDS")
	params:
		pheno  = config["phenotype_files"],
		covars = config["covariates_file"],
		script = "workflow/scripts/03-2_haplotypes_phased.R",
		haplo = lambda wc, input: input.prefix.replace(".sentinel", ".dosage.tsv"),
	log:
		ws_path("logs/build_haplotypes/{locus}_{dataset}_phased.log")
	conda:
		"../envs/environment.yml"
	benchmark:
		ws_path("benchmarks/execute_glm/{locus}_{dataset}.tsv")
	resources:
		runtime=lambda wc, attempt: attempt * 6000,
	shell:
		"""
		Rscript {params.script}  \
			--haplo {params.haplo}  \
            --pheno {params.pheno}  \
			--covariate {params.covars} \
			--output {output.result} 2> {log}
		"""
