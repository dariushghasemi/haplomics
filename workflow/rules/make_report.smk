
rule render_report:
	input:
		markdown = "workflow/scripts/04-0_report.qmd",
		plt_heat = lambda wc: ws_path(f"plot_heatmaps/{wc.locus}_{wc.dataset}_plot_heatmap.png"),
		plt_hapl = lambda wc: ws_path(f"plot_haplotypes/{wc.locus}_{wc.dataset}_plot_haplotypes.png"),
		res_rds  = lambda wc: ws_path(f"results/{wc.locus}_{wc.dataset}_tidied.RDS"),
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
		res_tidy = lambda wc: full_path(f"result_associations/{wc.locus}_{wc.dataset}_tidied.RDS"),
		hap_char = lambda wc: full_path(f"report/{wc.locus}_{wc.dataset}_characteristics.tsv"),
		hap_log  = lambda wc: full_path(f"logs/plot_haplotypes/{wc.locus}_{wc.dataset}.log"),
	conda:
		"../envs/environment.yml"
	threads: 1
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		mkdir -p {params.odir}
		
        # copy the shared .qmd to a temporary .qmd per job
        TEMP_QMD="{params.odir}/{params.locus}_{params.assay}.qmd"
        cp {input.markdown} $TEMP_QMD
		
		echo "Rendering report for LOCUS={params.locus} ASSAY={params.assay}"

		Rscript -e '
		rmarkdown::render(
		    input = "{params.odir}/{params.locus}_{params.assay}.qmd",
		    output_file = "{params.html}",
			output_dir  = "{params.odir}",
			clean  = TRUE,
		    params = list(
				LOCUS = "{params.locus}",
				ASSAY = "{params.assay}",
				hist = "{params.hist_abs}",
				anot = "{params.anot_abs}",
				hap1 = "{params.hap1_abs}",
				hap2 = "{params.hap2_abs}",
				heat = "{params.heat_abs}",
				summ = "{params.tbl_summ}",
				res  = "{params.res_tidy}",
				char = "{params.hap_char}",
				log  = "{params.hap_log}"
			)
		)
		'
		rm -f $TEMP_QMD
		"""