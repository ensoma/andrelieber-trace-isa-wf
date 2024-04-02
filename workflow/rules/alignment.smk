###################
## BWA Alignment ##
###################

# Index the combined payload-genome with BWA.
rule bwa_index_combined_payload_and_genome_fasta:
    input:
        fasta = outdir / 'genome/combined_payload_and_genome/{payload_genome}.fa'
    output:
        directory(outdir / 'genome/bwa_index_payload_and_genome/{payload_genome}')
    params:
        block_size = 25000000
    log: outdir / 'logs/{payload_genome}.bwa_index_combined_payload_and_genome_fasta.log'
    benchmark: outdir / 'benchmarks/{payload_genome}.bwa_index_combined_payload_and_genome_fasta.benchmark.txt'
    conda: '../envs/bwa.yaml'
    container: 'docker://quay.io/biocontainers/bwa:0.7.17--he4a0461_11'
    shell:
        """
        (
        mkdir -p {output}

        bwa index \
            -p {output}/index \
            -b {params.block_size} \
            {input}
        ) 2> {log}
        """

# Align the reads to the combined genome with bwa mem.
rule bwa_mem_aligment_to_payload_and_genome_index:
    input:
        R1 = outdir / 'trimmed/{sample_name}_R1.ir.readthrough.trimmed.fastq.gz',
        R2 = outdir / 'trimmed/{sample_name}_R2.ir.readthrough.trimmed.fastq.gz',
        index = lambda wildcards: outdir / f"genome/bwa_index_payload_and_genome/{genome_map[sample_map[wildcards.sample_name]]['payload_fasta_name']}_{genome_map[sample_map[wildcards.sample_name]]['genome_fasta_name']}"
    output:
        outdir / 'alignments/{sample_name}.aligned.sam'
    log: outdir / 'logs/{sample_name}.bwa_mem_aligment_to_payload_and_genome_index.log'
    benchmark: outdir / 'benchmarks/{sample_name}.bwa_mem_aligment_to_payload_and_genome_index.benchmark.txt'
    threads: 6
    conda: '../envs/bwa.yaml'
    container: 'docker://quay.io/biocontainers/bwa:0.7.17--he4a0461_11'
    shell:
        """
        (
        bwa mem \
            -t {threads} \
            -o {output} \
            {input.index}/index \
            {input.R1} \
            {input.R2}
        ) 2> {log}
        """

# Use samtools to get alignment stats.
rule stats_for_alignment_to_combined_payload_and_genome_index:
    input:
        [outdir / f"alignments/{sample_name}.aligned.sam" for sample_name in sample_map.keys()]
    output:
        flagstat = [outdir / f"alignments/{sample_name}.aligned.sam.flagstat" for sample_name in sample_map.keys()],
        idxstat = [outdir / f"alignments/{sample_name}.aligned.sam.idxstat" for sample_name in sample_map.keys()]
    params:
        outdir = outdir
    log: outdir / 'logs/stats_for_alignment_to_combined_payload_and_genome_index.log'
    benchmark: outdir / 'benchmarks/stats_for_alignment_to_combined_payload_and_genome_index.benchmark.txt'
    threads: 4
    conda: '../envs/samtools.yaml'
    container: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_1'
    shell:
        """
        (
        for SAM in {input}; do
            SAMPLE_NAME=$(basename ${{SAM}} .aligned.sam)

            samtools flagstat -@ {threads} ${{SAM}} \
                > {params.outdir}/alignments/${{SAMPLE_NAME}}.aligned.sam.flagstat

            samtools sort -@ {threads} ${{SAM}} | samtools idxstat -@ {threads} - \
                > {params.outdir}/alignments/${{SAMPLE_NAME}}.aligned.sam.idxstat
        done
        ) 2> {log}
        """