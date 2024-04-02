###########################
## Filter the Alignments ##
###########################

# Filter reads based on quality/flags, remove duplicates,
# and remove reads with too large of an R1 5' softclip.
# 1) Flag filtering:
#   - Remove unmapped reads
#   - Excluded flags: 3852 (secondary and supplementary alignments, read or mate unmapped, PCR/optical duplicates, paltform/vendor QC)
#   - Required flags: 3 (paired and properly paired)
#   - Map quality ≥ 30
#   - No alternatate hits reported by bwa mem (XA:Z:*) [Default: Off]
#   - No supplementary alignments (SA:Z:*) [Default: Off]
#   - Read is not orphaned after filtering
# 2) Remove PCR duplicates.
# 3) Remove reads with an R1 5' softclip > 5.
rule filter_and_position_sort_aligned_sam:
    input:
        outdir / 'alignments/{sample_name}.aligned.sam'
    output:
        # Flag filtered.
        bam_flag_filtered = outdir / 'alignments/{sample_name}.filtered.bam',
        bai_flag_filtered = outdir / 'alignments/{sample_name}.filtered.bam.bai',
        flagstat_flag_filtered = outdir / 'alignments/{sample_name}.filtered.bam.flagstat',
        idxstat_flag_filtered = outdir / 'alignments/{sample_name}.filtered.bam.idxstat',
        # Deduplicated.
        bam_deduped = outdir / 'alignments/{sample_name}.deduped.bam',
        bai_deduped = outdir / 'alignments/{sample_name}.deduped.bam.bai',
        flagstat_deduped = outdir / 'alignments/{sample_name}.deduped.bam.flagstat',
        idxstat_deduped = outdir / 'alignments/{sample_name}.deduped.bam.idxstat',
        # 5' softclip filtered.
        bam_softclip = outdir / 'alignments/{sample_name}.softclip_removed.bam',
        bai_softclip = outdir / 'alignments/{sample_name}.softclip_removed.bam.bai',
        flagstat_softclip = outdir / 'alignments/{sample_name}.softclip_removed.bam.flagstat',
        idxstat_softclip = outdir / 'alignments/{sample_name}.softclip_removed.bam.idxstat'
    params:
        max_memory = '768M',
        excluded_flags = 3852,
        required_flags = 3,
        min_mapq = 30,
        remove_if_any_alt = "grep -Ev '\s[XS]A:Z:\S+' |" if remove_if_any_alt_alignments else '',
        softclip_threshold = 5
    log: outdir / 'logs/{sample_name}.filter_and_position_sort_aligned_sam.log'
    benchmark: outdir / 'benchmarks/{sample_name}.filter_and_position_sort_aligned_sam.benchmark.txt'
    threads: 4
    conda: '../envs/is_softclip.yaml'
    container: 'docker://ensoma/is_softclip:1.4.0'
    shell:
        """
        (
        # Filter based on flags and quality scores.
        samtools view -@ {threads} -h -F {params.excluded_flags} -f {params.required_flags} -q {params.min_mapq} {input} |
            # Removes reads with reported alternate and supplementary alignments.
            {params.remove_if_any_alt}
            # To remove orphaned reads: name sort, fix mates, and then keep only paired reads.
            samtools sort -n -@ {threads} -m {params.max_memory} -u - |
            samtools fixmate -@ {threads} -m -u - - |
            samtools view -@ {threads} -f 2 -u - |
            # Position sort and remove duplicates.
            samtools sort \
                -@ {threads} -m {params.max_memory} -O BAM --write-index \
                -o {output.bam_flag_filtered}##idx##{output.bai_flag_filtered}  -

        # Flagstat and idxstat the flag filtered bam.
        samtools flagstat -@ {threads} {output.bam_flag_filtered} > {output.flagstat_flag_filtered}
        samtools idxstat -@ {threads} {output.bam_flag_filtered} > {output.idxstat_flag_filtered}

        # Remove PCR duplicates.
        samtools sort -n -@ {threads} -m {params.max_memory} -u {output.bam_flag_filtered} |
            samtools fixmate -@ {threads} -m -u - - |
            samtools sort -@ {threads} -m {params.max_memory} -u - |
            samtools markdup -@ {threads} -r -O BAM - {output.bam_deduped}
        
        # Index the deduplicated BAM file.
        samtools index -b -o {output.bai_deduped} {output.bam_deduped}

        # Return flag and idxstats for the deduplicated bam.
        samtools flagstat -@ {threads} {output.bam_deduped} > {output.flagstat_deduped}
        samtools idxstat -@ {threads} {output.bam_deduped} > {output.idxstat_deduped}

        # Remove too large of 5' R1 softclip.
        python3 workflow/scripts/is_softclip.py softclip_filter -i {output.bam_deduped} -s {params.softclip_threshold} |
            # Remove orphaned R2 reads.
            samtools sort -n -@ {threads} -m {params.max_memory} -u - |
            samtools fixmate -@ {threads} -u - - |
            samtools view -f3 -u -@ {threads} - |
            # Position sort the further filtered bam.
            samtools sort \
                -@ {threads} -O BAM --write-index -m {params.max_memory} \
                -o {output.bam_softclip}##idx##{output.bai_softclip} -
        
        # Flagstat and idxstat of the too-large 5' softclip filter.
        samtools flagstat -@ {threads} {output.bam_softclip} > {output.flagstat_softclip}
        samtools idxstat -@ {threads} {output.bam_softclip} > {output.idxstat_softclip}
        ) 2> {log}
        """