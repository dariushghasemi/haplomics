
rule get_locus;
	input:
	   script = "01_subsetting_vcf.sh",
	   position = "locus_window.txt"
	output:
	   file = "genotype/${locus}.vcf.gz"
	params: 
	   file = "${locus}.vcf.gz"
	shell:
	   """
	   {input.script}{params.file}
	   """
	 
rule get_dosage;
	input:
	   script = "03_extracting_dosage.sh"
	output:
	   file = "genotype/${locus}_dosage.txt"
	   
	params: 
	   file = ".gz"
	shell:
	   """
	   {input.script}{params.file}
	   """

