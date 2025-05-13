
rule get_locus:
	input:
		vcf = config["path_vcf"],
	output:
		ofile = "results/genotype/{locus}.vcf.gz",
	conda:
		"../envs/environment.yml"
	params:
		region = lambda wildcards: get_region(wildcards.locus),
	resources:
		runtime=lambda wc, attempt: attempt * 30,
	shell:
		"""
		echo "Region: {params.region}"
		
		# Check if the file exists and is readable
		if [ ! -f {input.vcf} ] || [ ! -r {input.vcf} ]; then
			echo "The file does not exist or cannot be read."
			exit 2
		# Check if the file is empty
		elif [ ! -s {input.vcf} ]; then
			echo "The file is empty."
			exit 3
		else
			echo "VCF must be compressed."
			bcftools view {input.vcf} -r {params.region} -Oz -o {output.ofile}
			bcftools index {output.ofile}
		fi
		"""