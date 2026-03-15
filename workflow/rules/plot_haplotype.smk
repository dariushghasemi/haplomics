
rule plot_haplotypes:
	input:
		script = "workflow/scripts/03-3_haplotypes_plot.R",
		rds    = rules.build_haplotypes.output.result,
		annot  = rules.annotate_variants.output.annotation,
		bim    = rules.get_dosage.output.bim,
	output:
		char = ws_path("report/{locus}_{dataset}_characteristics.tsv"),
		plt1 = ws_path("plot_haplotypes/{locus}_{dataset}_plot_haplotypes.png"),
		plt2 = ws_path("plot_haplotypes/{locus}_{dataset}_plot_haplotypes_shrinked.png")
	conda:
		"../envs/environment.yml"
	log:
		ws_path("logs/plot_haplotypes/{locus}_{dataset}.log")
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		Rscript {input.script}  \
			--rds      {input.rds}  \
			--annotation {input.annot}  \
			--variants {input.bim}  \
			--char     {output.char}  \
			--plotfull {output.plt1}  \
			--plottiny {output.plt2} 2> {log}
		"""