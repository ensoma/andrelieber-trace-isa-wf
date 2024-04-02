# Generate insertion site annotation plots.
# This is an R script that outputs 3 plots:
#   1) Genomic feature distribution of insertion sites.
#       - Promoter region is defined by the 'promoter' parameter in the rule 'annotate_insertion_sites'.
#       - By default the promoter is defined as -2000 to 200.
#       - Promoter region is also applicable to the promoter distribution plot.
#   2) Distribution of insertion sites around promoters.
#       - The number of bins for the histogram is defined by the parameter 'promoter_plot_bins'.
#       - It's set to 100 by default.
#   3) Top n insertion sites for each sample.
#       - top n sites as defined by the parameter 'n_top_sites'.
#       - Will be annotated with the gene name if the insertion site is genic.
# Additionally, a fourth plot will be generated that contains the insertion site fractions.
rule insertion_site_annotation_plots:
    input:
        [outdir / f"annotated_insertion_sites/{sample_name}.insertion_sites.annotated.tsv" for sample_name in sample_map.keys()]
    output:
        feature_dist_table = outdir / 'annotated_insertion_sites/is_feature_distribution_table.tsv',
        feature_dist_plot = outdir / 'annotated_insertion_sites/is_feature_distribution_plot.png',
        promoter_dist_plot = outdir / 'annotated_insertion_sites/promoter_distribution_plot.png',
        top_is_plot = outdir / 'annotated_insertion_sites/top_insertion_sites_plot.png',
        top_is_table = outdir / 'annotated_insertion_sites/top_insertion_sites_table.tsv',
        fracs_plot = outdir / 'insertion_sites/is_fractions_plot.png'
    params:
        outdir = outdir / 'annotated_insertion_sites',
        promoter_plot_bins = 100,
        n_top_sites = 10,
        dpi = 300
    log: outdir / 'logs/insertion_site_annotation_plots.log'
    benchmark: outdir / 'benchmarks/insertion_site_annotation_plots.benchmark.txt'
    conda: '../envs/is_annotations.yaml'
    container: 'docker://ensoma/is_anno:1.5.0'
    shell:
        """
        (
        Rscript workflow/scripts/is_annotations.R genomic_feature_plot \
            -i $(echo "{input}" | tr ' ' ',') \
            -o {output.feature_dist_plot} \
            -t {output.feature_dist_table} \
            --decreasing \
            -p {params.dpi}
        
        Rscript workflow/scripts/is_annotations.R promoter_plot \
            -i $(echo "{input}" | tr ' ' ',') \
            -o {output.promoter_dist_plot} \
            -b {params.promoter_plot_bins} \
            --decreasing \
            -p {params.dpi}
        
        Rscript workflow/scripts/is_annotations.R top_sites_plot \
            -i $(echo "{input}" | tr ' ' ',') \
            -o {output.top_is_plot} \
            -t {output.top_is_table} \
            -n {params.n_top_sites} \
            --decreasing \
            -p {params.dpi}
        
        Rscript workflow/scripts/is_annotations.R fractions \
            -i $(echo "{input}" | tr ' ' ',') \
            -o {output.fracs_plot} \
            -d {params.dpi}
        ) 2> {log}
        """