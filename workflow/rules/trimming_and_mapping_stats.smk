################################
## Trimming and Mapping Stats ##
################################

# Create a trimming and mapping stats file.
rule trimming_and_mapping_stats:
    input:
        unfiltered_R1 = fastq_dir / '{sample_name}_R1_001.fastq.gz',
        adapter_trimmed_R1 = outdir / 'trimmed/{sample_name}_R1.adapter.trimmed.fastq.gz',
        ir_trimmed_R1 = outdir / 'trimmed/{sample_name}_R1.ir.trimmed.fastq.gz',
        adapter_readthrough_trimmed_R1 = outdir / 'trimmed/{sample_name}_R1.adapter.readthrough.trimmed.fastq.gz',
        ir_readthrough_trimmed_R1 = outdir / 'trimmed/{sample_name}_R1.ir.readthrough.trimmed.fastq.gz',
        aligned_sam = outdir / 'alignments/{sample_name}.aligned.sam',
        filtered_bam = outdir / 'alignments/{sample_name}.filtered.bam',
        filtered_bai = outdir / 'alignments/{sample_name}.filtered.bam.bai',
        deduped_bam = outdir / 'alignments/{sample_name}.deduped.bam',
        deduped_bai = outdir / 'alignments/{sample_name}.deduped.bam.bai',
        softclip_removed_bam = outdir / 'alignments/{sample_name}.softclip_removed.bam',
        softclip_removed_bai = outdir / 'alignments/{sample_name}.softclip_removed.bam.bai',
        insertion_sites = outdir / 'insertion_sites/{sample_name}.insertion_sites.filtered.bed'
    output:
        outdir / 'stats/{sample_name}.trimming_and_mapping_stats.tsv'
    log: outdir / 'logs/{sample_name}.trimming_and_mapping_stats.log'
    benchmark: outdir / 'benchmarks/{sample_name}.trimming_and_mapping_stats.benchmark.txt'
    threads: 4
    conda: '../envs/samtools.yaml'
    container: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_1'
    shell:
        """
        (
        # File header.
        printf 'stat\tvalue\n' > {output}.tmp

        # Unfiltered R1 read count.
        printf "fastq_unfiltered\t%s\n" \
            $(gzip -dc {input.unfiltered_R1} | wc -l | awk '{{print $1/4}}') \
            >> {output}.tmp
        
        # Adapter trimmed R1 read count.
        printf "fastq_adapter_trimmed\t%s\n" \
            $(gzip -dc {input.adapter_trimmed_R1} | wc -l | awk '{{print $1/4}}') \
            >> {output}.tmp
        
        # IR trimmed R1 read count.
        printf "fastq_ir_trimmed\t%s\n" \
            $(gzip -dc {input.ir_trimmed_R1} | wc -l | awk '{{print $1/4}}') \
            >> {output}.tmp
        
        # Adapter readthrough trimmed R1 read count.
        printf "fastq_adapter_readthrough_trimmed\t%s\n" \
            $(gzip -dc {input.adapter_readthrough_trimmed_R1} | wc -l | awk '{{print $1/4}}') \
            >> {output}.tmp

        # IR readthrough trimmed R1 read count.
        printf "fastq_ir_readthrough_trimmed\t%s\n" \
            $(gzip -dc {input.ir_readthrough_trimmed_R1} | wc -l | awk '{{print $1/4}}') \
            >> {output}.tmp
        
        # Aligned read count.
        printf "sam_aligned\t%s\n" \
            $(samtools view -@ {threads} -F 12 {input.aligned_sam} | cut -f 1 -d $'\t' | sort -u | wc -l) \
            >> {output}.tmp
        
        # Filtered read count.
        printf "bam_filtered\t%s\n" \
            $(samtools view -@ {threads} -F 12 {input.filtered_bam} | cut -f 1 -d $'\t' | sort -u | wc -l) \
            >> {output}.tmp
        
        # Deduped read count.
        printf "bam_deduped\t%s\n" \
            $(samtools view -@ {threads} -F 12 {input.deduped_bam} | cut -f 1 -d $'\t' | sort -u | wc -l) \
            >> {output}.tmp
        
        # Softclip removed read count.
        printf "bam_softclip_removed\t%s\n" \
            $(samtools view -@ {threads} -F 12 {input.softclip_removed_bam} | cut -f 1 -d $'\t' | sort -u | wc -l) \
            >> {output}.tmp

        # Insertion site-associated reads after IS blacklist filtering.
        printf "blacklist_filtered_insertion_sites\t%s\n" \
            $(awk 'BEGIN{{OFS=FS="\t"; total=0}} {{total += $5}} END {{print total}}' {input.insertion_sites}) \
            >> {output}.tmp

        # Finish the output.
        awk \
            -v sample_name=$(echo {wildcards.sample_name} | sed -E 's/_S[0-9]+_L00[0-9]//') \
            'BEGIN{{FS=OFS="\t"}} NR == 1 {{print "sample_name",$0}} NR > 1 {{print sample_name,$0}}' \
            {output}.tmp \
            > {output}
        
        rm {output}.tmp
        ) 2> {log}
        """

# Aggregate the trimming and mapping stats.
rule aggregate_trimming_and_mapping_stats:
    input:
        [outdir / f"stats/{sample_name}.trimming_and_mapping_stats.tsv" for sample_name in sample_map.keys()]
    output:
        outdir / 'stats/aggregate.trimming_and_mapping_stats.tsv'
    log: outdir / 'logs/aggregate_trimming_and_mapping_stats.log'
    benchmark: outdir / 'benchmarks/aggregate_trimming_and_mapping_stats.benchmark.txt'
    conda: '../envs/csvtk.yaml'
    container: 'docker://quay.io/biocontainers/csvtk:0.28.0--h9ee0642_0'
    shell:
        """
        (
        csvtk concat -Tt {input} |
            csvtk spread -tT -k stat -v value |
            csvtk sort -Tt -k 1:N |
            csvtk mutate2 -Tt -w0 -n 'Adapter Trimmed (FASTQ)' -e '$fastq_unfiltered - $fastq_adapter_trimmed' |
            csvtk mutate2 -Tt -w0 -n 'IR Trimmed (FASTQ)' -e '$fastq_adapter_trimmed - $fastq_ir_trimmed' |
            csvtk mutate2 -Tt -w0 -n 'Adapter Readthrough Trimmed (FASTQ)' -e '$fastq_ir_trimmed - $fastq_adapter_readthrough_trimmed' |
            csvtk mutate2 -Tt -w0 -n 'IR Readthrough Trimmed (FASTQ)' -e '$fastq_adapter_readthrough_trimmed - $fastq_ir_readthrough_trimmed' |
            csvtk mutate2 -Tt -w0 -n 'Unmapped Reads (SAM)' -e '$fastq_ir_readthrough_trimmed - $sam_aligned' |
            csvtk mutate2 -Tt -w0 -n 'Flag Filtered (BAM)' -e '$sam_aligned - $bam_filtered' |
            csvtk mutate2 -Tt -w0 -n 'Deduplicated (BAM)' -e '$bam_filtered - $bam_deduped' |
            csvtk mutate2 -Tt -w0 -n "5' Softclip Filtered (BAM)" -e '$bam_deduped - $bam_softclip_removed' |
            csvtk mutate2 -Tt -w0 -n 'Blacklisted (ISs)' -e '$bam_softclip_removed - $blacklist_filtered_insertion_sites' |
            csvtk mutate2 -Tt -w0 -n 'Remaining' -e '$blacklist_filtered_insertion_sites' |
            csvtk cut -Tt -f -2--11 -o {output}
        ) 2> {log}
        """