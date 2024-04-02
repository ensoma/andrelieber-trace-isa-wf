###########################
## Payload Mapping Reads ##
###########################

# Extract the payload-mapping reads.
rule extract_payload_mapping_reads_from_alignment:
    input:
        bam = outdir / 'alignments/{sample_name}.softclip_removed.bam',
        bai = outdir / 'alignments/{sample_name}.softclip_removed.bam.bai',
        payload_fasta = lambda wildcards: outdir / f"genome/payload_only/{genome_map[sample_map[wildcards.sample_name]]['payload_fasta_name']}.fa"
    output:
        bam = outdir / 'alignments/{sample_name}.payload.bam',
        bai = outdir / 'alignments/{sample_name}.payload.bam.bai'
    log: outdir / 'logs/{sample_name}.extract_payload_mapping_reads_from_alignment.log'
    benchmark: outdir / 'benchmarks/{sample_name}.extract_payload_mapping_reads_from_alignment.benchmark.txt'
    threads: 4
    conda: '../envs/samtools.yaml'
    container: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_1'
    shell:
        """
        (
        samtools view \
            -@ {threads} \
            -O BAM \
            -o {output.bam}##idx##{output.bai} \
            --write-index \
            {input.bam} \
            $(head -n1 {input.payload_fasta} | sed 's/^>//')
        ) 2> {log}
        """