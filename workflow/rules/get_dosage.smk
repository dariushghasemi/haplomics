
rule get_dosage:
	input:
		script = "workflow/scripts/01-2_extracting_dosage.sh",
		vcf    = ws_path("genotype/{locus}.vcf.gz")
	output:
		dosage = ws_path("dosage/{locus}.dosage"),
		snps   = ws_path("dosage/{locus}.snps"),
		bim    = ws_path("dosage/{locus}.bim")
	conda:
		"../envs/environment.yml"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		# extract variants in the region after merging
		bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%DS\n]' {input.vcf} -o {output.dosage}
		bcftools query -f '%CHROM:%POS:%REF:%ALT\n' {input.vcf} -o {output.snps}
		bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' {input.vcf} -o {output.bim}
		"""
