
rule merge_data:
	input:
		script = "workflow/scripts/03-1_haplotypes_data.R",
		dosage = "results/dosage/{locus}.dosage",
		phenotype = lambda wildcards:  get_pheno(wildcards.dataset),
	output:
		odata = "results/merged_data/{locus}_{dataset}_merged_data.RDS"
		#summ  = "results/report/{locus}_merged_data_summary.txt"
	conda:
		"../envs/environment.yml"
	params:
		covariate = config.get("covariates_file"),
		min_ac = config.get("thresholds").get("min_ac"),
	log:
		"logs/merged_data/{locus}_{dataset}_merged_data.log"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		Rscript {input.script} \
			--dosage {input.dosage} \
			--phenotype {input.phenotype}  \
			--covariate {params.covariate}  \
			--min_ac {params.min_ac} \
			--output {output.odata}  2> {log}
		"""
