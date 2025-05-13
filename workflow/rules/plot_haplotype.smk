
rule plot_haplotypes:
	input:
		script = "workflow/scripts/03-3_haplotypes_plot.R",
		result = "results/result_associations/{locus}_{dataset}_association_results.RDS",
		annotation = "results/annotation/{locus}_summary.tsv",
		variants   = "results/dosage/{locus}.bim"
	output:
		plt1 = "results/plot_haplotypes/{locus}_{dataset}_plot_haplotypes.png",
		plt2 = "results/plot_haplotypes/{locus}_{dataset}_plot_haplotypes_shrinked.png"
	conda:
		"../envs/environment.yml"
	params:
		region="{locus}"
	log:
		"logs/plot_haplotypes/{locus}_{dataset}.log"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		Rscript {input.script}  \
			--rds {input.result}  \
			--annotation {input.annotation}  \
			--variants {input.variants}  \
			--locus {params.region}  \
			--output1 {output.plt1}  \
			--output2 {output.plt2} 2> {log}
		"""