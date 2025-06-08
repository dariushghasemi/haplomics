
rule render_report:
	input:
		markdown = "workflow/scripts/04-0_report.Rmd",
		plt_hist = "results/plot_histogram/{locus}.png",
		plt_anot = "results/plot_annotation/{locus}.png",
		plt_hap1 = lambda wc: expand("results/plot_haplotypes/{locus}_{dataset}_plot_haplotypes.png", locus = wc.locus, dataset = wc.dataset),
		plt_hap2 = lambda wc: expand("results/plot_haplotypes/{locus}_{dataset}_plot_haplotypes_shrinked.png",  locus = wc.locus, dataset = wc.dataset),
		plt_heat = lambda wc: expand("results/plot_heatmaps/{locus}_{dataset}_plot_heatmap.png", locus = wc.locus, dataset = wc.dataset),
		res_rds  = lambda wc: expand("results/result_tidied/{locus}_{dataset}_association_results_tidied.RDS", locus = wc.locus, dataset = wc.dataset),
		#tbl_summ = "results/report/{locus}_merged_data_summary.txt",
	output:
		html = "results/report_html/{locus}_{dataset}.nb.html"
	params:
		locus = "{locus}",
		assay = "{dataset}",
		html = "{locus}_{dataset}.nb.html",
		odir = "results/report_html/"
	conda:
		"../envs/environment.yml"
	threads: 1
	resources:
		runtime=lambda wc, attempt: attempt * 30,
		mem_mb=4048, disk_mb=20000
	shell:
		"""
		mkdir -p {params.odir}
		echo 'Start to create HTML report...'
		
		Rscript -e '
		rmarkdown::render(
		    input = "{input.markdown}",
		    output_file = "{params.html}",
			output_dir  = "{params.odir}",
		    params = list(
				LOCUS = "{params.locus}",
				ASSAY = "{params.assay}",
				hist = "{input.plt_hist}",
				anot = "{input.plt_anot}",
				hap1 = "{input.plt_hap1}",
				hap2 = "{input.plt_hap2}",
				heat = "{input.plt_heat}", 
				res  = "{input.res_rds}"
			)
		)
		'
		"""