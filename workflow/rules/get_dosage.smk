
rule get_dosage:
	input:
		script = "workflow/scripts/01-2_extracting_dosage.sh",
		vcf    = "results/genotype/{locus}.vcf.gz"
	output:
		dosage = "results/dosage/{locus}.dosage",
		snps   = "results/annotation/{locus}.snps",
		#annotation = "results/annotation/{locus}_annotation.txt"
	conda:
		"../envs/environment.yml"
	resources:
		runtime=lambda wc, attempt: attempt * 30,
		mem_mb=get_mem_plt, disk_mb=5000
	shell:
		"""
		# extract variants in the region after merging
		bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%DS\n]' {input.vcf} -o {output.dosage}
		bcftools query -f '%CHROM:%POS:%REF:%ALT\n' {input.vcf} -o {output.snps}
		"""
