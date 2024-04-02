####################
## Prepare Genome ##
####################

# harmonize payload fasta files.
#   - If the file is gzipped, unzip it.
#   - If the file is not gzipped, move it.
#   - Index the fasta file.
rule harmonize_and_index_payload_fasta:
    input:
        lambda wildcards: genome_dir / f"{payload_name_file[wildcards.payload_name]}"
    output:
        fasta = outdir / 'genome/payload_only/{payload_name}.fa',
        fai = outdir / 'genome/payload_only/{payload_name}.fa.fai'
    log: outdir / 'logs/{payload_name}.harmonize_and_index_payload_fasta.log'
    benchmark: outdir / 'benchmarks/{payload_name}.harmonize_and_index_payload_fasta.benchmark.txt'
    conda: '../envs/samtools.yaml'
    container: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_1'
    shell:
        """
        (
        if [[ {input} =~ \.gz$ ]]; then
            gzip -dc {input} > {output.fasta}
        else
            cp {input} {output.fasta}
        fi

        samtools faidx -o {output.fai} {output.fasta}
        ) 2> {log}
        """

# harmonize genome fasta files.
#   - If the file is gzipped, unzip it.
#   - If the file is not gzipped, move it.
#   - Index the fasta file.
rule harmonize_and_index_genome_fasta:
    input:
        lambda wildcards: genome_dir / f"{genome_name_file[wildcards.genome_name]}"
    output:
        fasta = outdir / 'genome/genome_only/{genome_name}.fa',
        fai = outdir / 'genome/genome_only/{genome_name}.fa.fai'
    log: outdir / 'logs/{genome_name}.harmonize_and_index_genome_fasta.log'
    benchmark: outdir / 'benchmarks/{genome_name}.harmonize_and_index_genome_fasta.benchmark.txt'
    conda: '../envs/samtools.yaml'
    container: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_1'
    shell:
        """
        (
        if [[ {input} =~ \.gz$ ]]; then
            gzip -dc {input} > {output.fasta}
        else
            cp {input} {output.fasta}
        fi

        samtools faidx -o {output.fai} {output.fasta}
        ) 2> {log}
        """

# Combine the payload and genome and index it.
rule combine_payload_and_genome_fastas:
    input:
        payload_fasta = lambda wildcards: outdir / f"genome/payload_only/{payload_genome_combos[wildcards.payload_genome]['payload_fasta_name']}.fa",
        genome_fasta = lambda wildcards: outdir / f"genome/genome_only/{payload_genome_combos[wildcards.payload_genome]['genome_fasta_name']}.fa"
    output:
        fasta = outdir / 'genome/combined_payload_and_genome/{payload_genome}.fa',
        fai = outdir / 'genome/combined_payload_and_genome/{payload_genome}.fa.fai'
    log: outdir / 'logs/{payload_genome}.combine_payload_and_genome_fastas.log'
    benchmark: outdir / 'benchmarks/{payload_genome}.combine_payload_and_genome_fastas.benchmark.txt'
    conda: '../envs/samtools.yaml'
    container: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_1'
    shell:
        """
        (
        cat {input.genome_fasta} {input.payload_fasta} > {output.fasta}
        samtools faidx -o {output.fai} {output.fasta}
        ) 2> {log}
        """