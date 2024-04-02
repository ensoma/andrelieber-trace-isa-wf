#################
# Sequence Logo #
#################

# Expand the IS ranges on each side and retrieve the sequences.
# Distance is set by the 'slop_size' parameter (20 by default).
rule retrieve_insertion_site_sequences:
    input:
        bed = outdir / 'insertion_sites/{sample_name}.insertion_sites.filtered.bed',
        genome_fasta = lambda wildcards: outdir / f"genome/combined_payload_and_genome/{genome_map[sample_map[wildcards.sample_name]]['payload_fasta_name']}_{genome_map[sample_map[wildcards.sample_name]]['genome_fasta_name']}.fa",
        genome_fai = lambda wildcards: outdir / f"genome/combined_payload_and_genome/{genome_map[sample_map[wildcards.sample_name]]['payload_fasta_name']}_{genome_map[sample_map[wildcards.sample_name]]['genome_fasta_name']}.fa.fai"
    output:
        outdir / 'seqlogo/{sample_name}.insertion_site_seq.fasta'
    params:
        slop_size = 20
    log: outdir / 'logs/{sample_name}.retrieve_insertion_site_sequences.log'
    benchmark: outdir / 'benchmarks/{sample_name}.retrieve_insertion_site_sequences.benchmark.txt'
    conda: '../envs/bedtools.yaml'
    container: 'docker://quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_3'
    shell:
        """
        (
        bedtools slop -b {params.slop_size} -g {input.genome_fai} -i {input.bed} |
            bedtools getfasta -s -fi {input.genome_fasta} -fo {output} -bed -
        ) 2> {log}
        """

# Generate the sequence logo.
# Uses the R ggseqlogo library.
rule create_sequence_logo:
    input:
        [outdir / f"seqlogo/{sample_name}.insertion_site_seq.fasta" for sample_name in sample_map.keys()]
    output:
        outdir / 'seqlogo/is_seqlogo.png'
    params:
        dpi = 300,
        height_multiplier = 1.25
    log: outdir / 'logs/create_sequence_logo.log'
    benchmark: outdir / 'benchmarks/create_sequence_logo.benchmark.txt'
    conda: '../envs/is_annotations.yaml'
    container: 'docker://ensoma/is_anno:1.5.0'
    shell:
        """
        (
        Rscript workflow/scripts/is_annotations.R seqlogo \
            --input=$(echo "{input}" | tr ' ' ',') \
            --output={output} \
            --dpi={params.dpi} \
            --height_multiplier={params.height_multiplier}
        ) 2> {log}
        """