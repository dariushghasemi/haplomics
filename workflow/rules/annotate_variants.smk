
rule annotate_variants:
	input:
		vep_api = "workflow/scripts/01-6_vep_annotator.py",
		snps    = ws_path("dosage/{locus}.snps"),
	output:
		annotation = ws_path("annotation/{locus}_summary.tsv"),
		json       = ws_path("annotation/{locus}_output.json")
	conda:
		"../envs/environment.yml"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		source /exchange/healthds/singularity_functions
		module load python-3.9.10/py-requests/2.31.0
		module load python-3.9.10/py-pandas/1.5.3

		# run the VEP annotator
		python {input.vep_api} --snps-list {input.snps} --output-summary {output.annotation} --output-json {output.json}
		"""
