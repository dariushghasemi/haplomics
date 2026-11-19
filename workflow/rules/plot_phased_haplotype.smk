
rule plot_haplotypes:
    input:
        freqs = rules.compute_matrix.output.frequencies,
        variants = rules.compute_matrix.output.variants,
        script = "workflow/scripts/haplotypes_plot_phased.py"
    output:
        pdf = ws_path("plots/{locus}.haplotypes.pdf"),
        png = ws_path("plots/{locus}.haplotypes.png"),
        svg = ws_path("plots/{locus}.haplotypes.svg"),
        manifest = ws_path("plots/{locus}.haplotypes.manifest.json"),
        sentinel = touch(ws_path("plots/{locus}.haplotypes.sentinel")),
    params:
        prefix = lambda wc, output: output.sentinel.replace(".sentinel", ""),
        max_haplotypes=30,
        min_frequency=0.0,
        rows_per_page=30,
        variants_per_page=30
    conda:
        "../envs/python.yaml"
    #log:
    #    ws_path("logs/plot_haplotypes/{locus}.log")
    benchmark:
        ws_path("benchmarks/plot_haplotype/{locus}.tsv")
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    threads: 1
    shell:
        """
		python {input.script:q} \
          --frequencies {input.freqs:q} \
		  --variants {input.variants:q} \
          --output-prefix {params.prefix:q} \
          --max-haplotypes {params.max_haplotypes} \
          --min-frequency {params.min_frequency} \
          --rows-per-page {params.rows_per_page} \
          --variants-per-page {params.variants_per_page}
        """
#  2> {log}