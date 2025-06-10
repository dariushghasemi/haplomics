import os

rule render_report:
	input:
		markdown = "workflow/scripts/04-0_report.qmd",
		plt_heat = ws_path("plot_heatmaps/{locus}_{dataset}_plot_heatmap.png"),
		res_rds  = ws_path("result_tidied/{locus}_{dataset}_association_results_tidied.RDS"),
	output:
		html = ws_path("report_html/{locus}_{dataset}.nb.html")
	params:
		locus = "{locus}",
		assay = "{dataset}",
		html = "{locus}_{dataset}.nb.html",
		odir = ws_path("report_html/"),
		hist_abs = lambda wc: full_path(f"plot_histogram/{wc.locus}.png"),
		anot_abs = lambda wc: full_path(f"plot_annotation/{wc.locus}.png"),
		hap1_abs = lambda wc: full_path(f"plot_haplotypes/{wc.locus}_{wc.dataset}_plot_haplotypes.png"),
		hap2_abs = lambda wc: full_path(f"plot_haplotypes/{wc.locus}_{wc.dataset}_plot_haplotypes_shrinked.png"),
		heat_abs = lambda wc: full_path(f"plot_heatmaps/{wc.locus}_{wc.dataset}_plot_heatmap.png"),
		tbl_summ = lambda wc: full_path(f"report/{wc.locus}_{wc.dataset}_merged_data_summary.txt"),
		res_abs  = lambda wc: full_path(f"result_tidied/{wc.locus}_{wc.dataset}_association_results_tidied.RDS")
	conda:
		"../envs/environment.yml"
	threads: 1
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		mkdir -p {params.odir}
		
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
				heat = "{params.heat_abs}",
				summ = "{params.tbl_summ}",
				res  = "{params.res_abs}"
			)
		)
		'
		"""