
rule plot_haplotypes:
	input:
		script = "workflow/scripts/03-3_haplotypes_plot.R",
		result = "results/result_associations/{locus}_{dataset}_association_results.RDS",
		#annotation = "results/annotation/{locus}_annotation.txt",
		#variants = "results/annotation/{locus}_variants.list"
	output:
		plt1 = "results/plot_haplotypes/{locus}_{dataset}_plot_haplotypes.png",
		plt2 = "results/plot_haplotypes/{locus}_{dataset}_plot_haplotypes_shrinked.png"
	conda:
		"../envs/environment.yml"
	params:
		#result
	log:
		"logs/plot_haplotypes/{locus}_{dataset}.log"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
		mem_mb=get_mem_plt, disk_mb=20000
	shell:
		"""
		Rscript {input.script} {input.result} {output.plt1} {output.plt2} 2> {log}
		"""
#{input.annotation} {input.variants} 