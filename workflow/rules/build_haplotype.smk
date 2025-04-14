
rule build_haplotypes:
	input:
		script = "workflow/scripts/03-2_haplotypes_building.R",
		data   = "results/merged_data/{locus}_{dataset}_merged_data.RDS"
	output:
		result = "results/result_associations/{locus}_{dataset}_association_results.RDS"
	conda:
		"../envs/environment.yml"
	params:
		#data
	log:
		"logs/build_haplotypes/{locus}_{dataset}.log"
	resources:
		runtime=lambda wc, attempt: attempt * 6000,
		mem_mb=get_mem_mb
	shell:
		"""
		Rscript {input.script} {input.data} {output.result} 2> {log}
		"""
