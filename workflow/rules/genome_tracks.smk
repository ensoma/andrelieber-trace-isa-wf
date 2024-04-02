####################
## Genomic Tracks ##
####################

# Convert the BAM files to bigwigs.
rule bam_to_bedgraph:
    input:
        bams = [outdir / f"alignments/{sample_name}.softclip_removed.bam" for sample_name in sample_map.keys()],
        bais = [outdir / f"alignments/{sample_name}.softclip_removed.bam.bai" for sample_name in sample_map.keys()]
    output:
        [outdir / f"tracks/{sample_name}.bedgraph" for sample_name in sample_map.keys()]
    params:
        outdir = outdir
    log: outdir / 'logs/bam_to_bedgraph.log'
    benchmark: outdir / 'benchmarks/bam_to_bedgraph.benchmark.txt'
    conda: '../envs/bedtools.yaml'
    container: 'docker://quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_3'
    shell:
        """
        (
        for BAM in {input.bams}; do
            SAMPLE_NAME=$(basename ${{BAM}} .softclip_removed.bam)

            bedtools genomecov -ibam ${{BAM}} -bg -pc |
                sort -k1,1 -k2,2n > {params.outdir}/tracks/${{SAMPLE_NAME}}.bedgraph
        done
        ) 2> {log}
        """

# Create bigwig files from the bedgraph files.
rule bedgraph_to_bigwig:
    input:
        bedgraph = outdir / 'tracks/{sample_name}.bedgraph',
        genome_fai = (
            lambda wildcards: outdir 
            / f"genome/combined_payload_and_genome/"
            f"{genome_map[sample_map[wildcards.sample_name]]['payload_fasta_name']}_"
            f"{genome_map[sample_map[wildcards.sample_name]]['genome_fasta_name']}.fa.fai"
        )
    output:
        outdir / 'tracks/{sample_name}.bw'
    log: outdir / 'logs/{sample_name}.bedgraph_to_bigwig.log'
    benchmark: outdir / 'benchmarks/{sample_name}.bedgraph_to_bigwig.benchmark.txt'
    conda: '../envs/bedgraphtobigwig.yaml'
    container: 'docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0'
    shell:
        """
        (
        # If the input bedgraph is empty create an empty bigwig file.
        # else convert it to a bigwig.
        if [[ $(wc -l < {input.bedgraph}) -eq 0 ]]; then
            touch {output}
        else
            bedGraphToBigWig \
                {input.bedgraph} \
                {input.genome_fai} \
                {output}
        fi
        ) 2> {log}
        """