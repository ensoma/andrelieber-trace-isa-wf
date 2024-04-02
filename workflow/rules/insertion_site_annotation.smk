###############################
## Insertion Site Annotation ##
###############################

# Annotate the insertion sites.
# Uses the R ChIPseeker library to annotate insertion sites relative to genomic features.
#   - The promoter region is defined by the 'promoter' parameter (-2000 to +200 by default).
#   - 'level' parameter controls whether annotations are relative to genes or transcripts (transcripts by defult).
rule annotate_insertion_sites:
    input:
        bed = outdir / 'insertion_sites/{sample_name}.insertion_sites.filtered.bed',
        gtf = lambda wildcards: genome_dir / f"{genome_map[sample_map[wildcards.sample_name]]['genome_anno_file']}"
    output:
        outdir / 'annotated_insertion_sites/{sample_name}.insertion_sites.annotated.tsv'
    params:
        promoter = '-2000,200',
        level = 'transcript'
    log: outdir / 'logs/{sample_name}.annotate_insertion_sites.log'
    benchmark: outdir / 'benchmarks/{sample_name}.annotate_insertion_sites.benchmark.txt'
    conda: '../envs/is_annotations.yaml'
    container: 'docker://ensoma/is_anno:1.5.0'
    shell:
        """
        (
        Rscript workflow/scripts/is_annotations.R annotate \
            --input={input.bed} \
            --genome={input.gtf} \
            --promoter='{params.promoter}' \
            --level={params.level} \
            --output={output}
        ) 2> {log}
        """

# Aggregate the insertion site annotations.
# Concatenate the annotated insertion sites and then sort by sample name and then insertion site frequency.
rule aggregate_insertion_site_annotations:
    input:
        [outdir / f"annotated_insertion_sites/{sample_name}.insertion_sites.annotated.tsv" for sample_name in sample_map.keys()]
    output:
        outdir / 'annotated_insertion_sites/aggregated_insertion_site_annotations.tsv'
    log: outdir / 'logs/aggregate_insertion_site_annotations.log'
    benchmark: outdir / 'benchmarks/aggregate_insertion_site_annotations.benchmark.txt'
    conda: '../envs/csvtk.yaml'
    container: 'docker://quay.io/biocontainers/csvtk:0.28.0--h9ee0642_0'
    shell:
        """
        (
        csvtk concat -Tt {input} |
            csvtk sort -k 1:N -k 7:nr -Tt -o {output}
        ) 2> {log}
        """