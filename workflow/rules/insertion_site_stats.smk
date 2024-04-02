############################
## Insertion Site Summary ##
############################

# Create a discarded IS report for multiqc.
rule discarded_is_report:
    input:
        [outdir / f"insertion_sites/{sample}.insertion_sites.discarded.bed" for sample in sample_map.keys()]
    output:
        outdir / 'insertion_sites/discarded_is_report.tsv'
    log: outdir / 'logs/discarded_is_report.log'
    benchmark: outdir / 'benchmarks/discarded_is_report.benchmark.txt'
    conda: '../envs/csvtk.yaml'
    container: 'docker://quay.io/biocontainers/csvtk:0.28.0--h9ee0642_0'
    shell:
        """
        (
        printf 'sample_name\tseqname:start-end:strand\tscore\n' > {output}.tmp
        sample_names=()
        for file in {input}; do
            sample_name=$(basename ${{file}} .insertion_sites.discarded.bed | sed -E 's/_S[0-9]+_L00[0-9]//')
            sample_names+=(${{sample_name}})

            awk \
                -v sample_name=${{sample_name}} \
                'BEGIN{{FS=OFS="\t"}} NR > 1 {{print sample_name,$1":"$2"-"$3":"$6,$5}}' \
                ${{file}} \
                >> {output}.tmp
        done

        # Pivot the report to wider if there is one or more samples.
        # If thare are not discarded insertion sites print the header.
        if [[ $(wc -l < {output}.tmp) -ge 2 ]]; then
            csvtk sort -k 2:N -Tt {output}.tmp |
                csvtk spread -k1 -v3 -Tt -o {output}
        else
            printf \
                'seqname:start-end:strand\t%s\n' "$(echo ${{sample_names[@]}} | tr ' ' '\t')" \
                > {output}
        fi

        rm {output}.tmp
        ) 2> {log}
        """

# Summary file of useful stats surrounding insertion sites.
rule insertion_site_stats:
    input:
        payload_bam = outdir / 'alignments/{sample_name}.payload.bam',
        payload_bai = outdir / 'alignments/{sample_name}.payload.bam.bai',
        insertion_sites = outdir / 'insertion_sites/{sample_name}.insertion_sites.filtered.bed',
        unfiltered_sites = outdir / 'insertion_sites/{sample_name}.insertion_sites.bed',
        blacklisted_sites = outdir / 'insertion_sites/{sample_name}.insertion_sites.discarded.bed',
        genome_fai = lambda wildcards: outdir / f"genome/genome_only/{genome_map[sample_map[wildcards.sample_name]]['genome_fasta_name']}.fa.fai"
    output:
        outdir / 'stats/{sample_name}.insertion_site_stats.tsv'
    params:
        payload_file_name = lambda wildcards: genome_map[sample_map[wildcards.sample_name]]['payload_fasta_name'],
        payload_start_loc = lambda wildcards: regions_map[genome_map[sample_map[wildcards.sample_name]]['payload_fasta_name']]['start_loc'],
        payload_end_loc = lambda wildcards: regions_map[genome_map[sample_map[wildcards.sample_name]]['payload_fasta_name']]['end_loc']
    log: outdir / 'logs/{sample_name}.insertion_site_stats.log'
    benchmark: outdir / 'benchmarks/{sample_name}.insertion_site_stats.benchmark.txt'
    conda: '../envs/samtools.yaml'
    container: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_1'
    shell:
        """
        (
        # Payload name in BAM file.
        payload_name=$(samtools view -@ {threads} {input.payload_bam} | cut -f3 -d$'\t' | uniq)

        # File header.
        printf 'stat\tvalue\n' > {output}.tmp

        # Insertion sites.
        insertion_sites=$(wc -l < {input.insertion_sites})
        printf "insertion_sites\t%s\n" ${{insertion_sites}} >> {output}.tmp

        # Blacklisted insertion sites.
        blacklisted_sites=$(wc -l < {input.blacklisted_sites})
        printf "blacklisted_sites\t%s\n" ${{blacklisted_sites}} >> {output}.tmp

        # Total read pairs.
        total_reads=$(awk -F'\t' 'BEGIN{{total=0}}{{total += $5}}END{{print total}}' {input.insertion_sites})
        printf "total_reads\t%s\n" ${{total_reads}} >> {output}.tmp

        # Blacklisted read pairs.
        blacklisted_reads=$(awk -F'\t' 'BEGIN{{total=0}}{{total += $5}}END{{print total}}' {input.blacklisted_sites})
        printf "blacklisted_reads\t%s\n" ${{blacklisted_reads}} >> {output}.tmp

        # Genomic read pairs.
        genomic_reads=$(awk -v payload_name=${{payload_name}} 'BEGIN{{OFS=FS="\t"; total=0}} $1 != payload_name {{total += $5}} END {{print total}}' {input.insertion_sites})
        printf "genome_reads\t%s\n" ${{genomic_reads}} >> {output}.tmp

        # Payload insertion sites.
        payload_insertion_sites=$(awk -v payload_name=${{payload_name}} 'BEGIN{{OFS=FS="\t"}} $1 == payload_name' {input.unfiltered_sites} | wc -l)
        printf "payload_insertion_sites\t%s\n" ${{payload_insertion_sites}} >> {output}.tmp

        # Payload read pairs.
        payload_reads=$(awk -v payload_name=${{payload_name}} 'BEGIN{{OFS=FS="\t"; total=0}} $1 == payload_name {{total += $5}} END {{print total}}' {input.unfiltered_sites})
        printf "payload_reads\t%s\n" ${{payload_reads}} >> {output}.tmp

        # Payload region pairs.
        payload_region=$(samtools view -c -f 64 {input.payload_bam} ${{payload_name}}:{params.payload_start_loc}-{params.payload_end_loc})
        printf "payload_region\t%s\n" ${{payload_region}} >> {output}.tmp
        
        # Fill out the rest of the file.
        awk \
            -v payload_name={params.payload_file_name} \
            -v sample_name=$(echo {wildcards.sample_name} | sed -E 's/_S[0-9]+_L00[0-9]//') \
            'BEGIN {{FS=OFS="\t"}} NR == 1 {{print "sample_name","payload_name",$0}} NR > 1 {{print sample_name,payload_name,$0}}' \
            {output}.tmp \
            > {output}

        rm {output}.tmp
        ) 2> {log}
        """

# Aggregate the stats.
rule aggregate_insertion_site_stats:
    input:
        insertion_site_stats = [
            outdir / f"stats/{sample_name}.insertion_site_stats.tsv" for sample_name in sample_map.keys()
        ],
    output:
        outdir / 'stats/aggregate_stats.tsv'
    log: outdir / 'logs/aggregate_stats.log'
    benchmark: outdir / 'benchmarks/aggregate_stats.benchmark.txt'
    conda: '../envs/csvtk.yaml'
    container: 'docker://quay.io/biocontainers/csvtk:0.28.0--h9ee0642_0'
    shell:
        """
        (
        # Concatenate the aggregate stats and pivot wider.
        csvtk concat -Tt {input} |
            csvtk spread -k stat -v value -tT |
            csvtk sort -k 1:N -tT -o {output}
        ) 2> {log}
        """