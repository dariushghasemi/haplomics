import os

rule render_report:
	input:
		markdown = "workflow/scripts/04-0_report.qmd",
		plt_hist = ws_path("plot_histogram/{locus}.png"),
		plt_anot = ws_path("plot_annotation/{locus}.png"),
		plt_hap1 = lambda wc: expand(ws_path("plot_haplotypes/{locus}_{dataset}_plot_haplotypes.png"), locus = wc.locus, dataset = wc.dataset),
		plt_hap2 = lambda wc: expand(ws_path("plot_haplotypes/{locus}_{dataset}_plot_haplotypes_shrinked.png"),  locus = wc.locus, dataset = wc.dataset),
		plt_heat = lambda wc: expand(ws_path("plot_heatmaps/{locus}_{dataset}_plot_heatmap.png"), locus = wc.locus, dataset = wc.dataset),
		res_rds  = lambda wc: expand(ws_path("result_tidied/{locus}_{dataset}_association_results_tidied.RDS"), locus = wc.locus, dataset = wc.dataset),
		#tbl_summ = "results/report/{locus}_merged_data_summary.txt",
	output:
		html = ws_path("report_html/{locus}_{dataset}.nb.html")
	params:
		locus = "{locus}",
		assay = "{dataset}",
		html = "{locus}_{dataset}.nb.html",
		odir = ws_path("report_html/"),
		hist_abs = lambda wc: os.path.abspath(f"results/plot_histogram/{wc.locus}.png"),
		anot_abs = lambda wc: os.path.abspath(f"results/plot_annotation/{wc.locus}.png"),
		hap1_abs = lambda wc: os.path.abspath(f"results/plot_haplotypes/{wc.locus}_{wc.dataset}_plot_haplotypes.png"),
		hap2_abs = lambda wc: os.path.abspath(f"results/plot_haplotypes/{wc.locus}_{wc.dataset}_plot_haplotypes_shrinked.png"),
		heat_abs = lambda wc: os.path.abspath(f"results/plot_heatmaps/{wc.locus}_{wc.dataset}_plot_heatmap.png"),
		res_abs  = lambda wc: os.path.abspath(f"results/result_tidied/{wc.locus}_{wc.dataset}_association_results_tidied.RDS")
	conda:
		"../envs/environment.yml"
	threads: 1
	resources:
		runtime=lambda wc, attempt: attempt * 30,
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
				hist = "{params.hist_abs}",
				anot = "{params.anot_abs}",
				hap1 = "{params.hap1_abs}",
				hap2 = "{params.hap2_abs}",
				heat = "{params.heat_abs}", 
				res  = "{params.res_abs}"
			)
		)
		'
		"""