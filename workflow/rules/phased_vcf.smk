
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


