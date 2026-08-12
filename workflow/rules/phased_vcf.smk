
rule get_phasing:
    input:
        rules.get_locus.output.ofile,
    output:
        gt = ws_path("phased/{locus}_phased_gt.tsv"),
        samples = ws_path("phased/{locus}.samples"),
    params:
    conda:
        "../envs/environment.yml"
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    shell:
        """
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' {input} > {output.gt}
        bcftools query -l {input} > {output.samples}
        """


rule compute_matrix:
    input:
        phasing = rules.get_phasing.output.gt,
        samples = rules.get_phasing.output.samples,
    output:
        sentinel = ws_path("phased/{locus}.sentinel"),
    params:
        script = "workflow/scripts/02-1_haplotype_matrix.py",
        min_freq = config.get("thresholds").get("min_freq"),
        prefix = lambda wc, output: output.sentinel.replace(".sentinel", ""),
    conda:
        "../envs/environment.yml"
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    shell:
        """
        source /exchange/healthds/singularity_functions
		module load python-3.9.10/py-requests/2.31.0
		module load python-3.9.10/py-pandas/1.5.3
        
        python {params.script} \
            --gt       {input.phasing} \
            --samples  {input.samples} \
            --min-freq {params.min_freq} \
            --out-prefix {params.prefix}
        
        touch {output.sentinel}
        """
