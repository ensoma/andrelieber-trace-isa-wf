##########################
# Insertion Site Calling #
##########################

# Create a BED file with IS scores.
# ISs within 5 bases are collapsed.
#   - The IS score is the sum of scores for positions within 5 bases.
#   - Position is the maximum score.
# Position sort the final bed.
rule count_insertion_sites:
    input:
        bams = [outdir / f"alignments/{sample_name}.softclip_removed.bam" for sample_name in sample_map.keys()],
        bais = [outdir / f"alignments/{sample_name}.softclip_removed.bam.bai" for sample_name in sample_map.keys()],
    output:
        [outdir / f"insertion_sites/{sample_name}.insertion_sites.bed" for sample_name in sample_map.keys()]
    params:
        outdir = outdir,
        grouping_distance = 5
    log: outdir / 'logs/count_insertion_sites.log'
    benchmark: outdir / 'benchmarks/count_insertion_sites.benchmark.txt'
    conda: '../envs/is_softclip.yaml'
    container: 'docker://ensoma/is_softclip:1.4.0'
    shell:
        """
        (
        for BAM in {input.bams}; do
            SAMPLE_NAME=$(basename ${{BAM}} .softclip_removed.bam)

            python3 workflow/scripts/is_softclip.py is_counter -i ${{BAM}} -g {params.grouping_distance} |
                sort -Vk1,2 \
                > {params.outdir}/insertion_sites/${{SAMPLE_NAME}}.insertion_sites.bed
        done
        ) 2> {log}
        """

# Create an insertion site blacklist BED file.
# Blacklist an insertion site if more than n samples have an IS score > 1.
rule insertion_site_blacklist_bed:
    input:
        [outdir / f"insertion_sites/{sample_name}.insertion_sites.bed" for sample_name in sample_map.keys()]
    output:
        outdir / 'insertion_sites/insertion_site_blacklist.bed'
    params:
        max_allowed_overlaps = 2,
        min_score = 2
    log: outdir / 'logs/insertion_site_blacklist_bed.log'
    benchmark: outdir / 'benchmarks/insertion_site_blacklist_bed.benchmark.txt'
    conda: '../envs/coreutils.yaml'
    container: 'docker://ubuntu:18.04'
    shell:
        """
        (
        cat {input} |
            awk -v min_score={params.min_score} 'BEGIN{{FS="\t"; OFS=" "}} $5 >= min_score {{print $1,$2,$3,$6}}' |
            sort |
            uniq -c |
            awk -v max={params.max_allowed_overlaps} 'BEGIN{{FS=" "; OFS="\t"}} $1 > max {{print $2,$3,$4,".",$1,$5}}' |
            sort -k1,1V -k2,3n -k6,6V \
            > {output}
        ) 2> {log}
        """

# Filter the insertion sites using the blacklist.
rule filter_blacklisted_insertion_sites:
    input:
        insertion_sites = [outdir / f"insertion_sites/{sample_name}.insertion_sites.bed" for sample_name in sample_map.keys()],
        blacklist = outdir / 'insertion_sites/insertion_site_blacklist.bed'
    output:
        retained = [outdir / f"insertion_sites/{sample_name}.insertion_sites.filtered.bed" for sample_name in sample_map.keys()],
        discarded = [outdir /f"insertion_sites/{sample_name}.insertion_sites.discarded.bed" for sample_name in sample_map.keys()]
    params:
        outdir = outdir
    log: outdir / 'logs/filter_blacklisted_insertion_sites.log'
    benchmark: outdir / 'benchmarks/filter_blacklisted_insertion_sites.benchmark.txt'
    conda: '../envs/bedtools.yaml'
    container: 'docker://quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_3'
    shell:
        """
        (
        for BED in {input.insertion_sites}; do
            SAMPLE_NAME=$(basename ${{BED}} .insertion_sites.bed)

            # Retained insertion sites.
            bedtools intersect \
                -s -v \
                -a ${{BED}} \
                -b {input.blacklist} \
                > {params.outdir}/insertion_sites/${{SAMPLE_NAME}}.insertion_sites.filtered.bed

            # Discarded insertion sites.
            bedtools intersect \
                -s -wa \
                -a ${{BED}} \
                -b {input.blacklist} \
                > {params.outdir}/insertion_sites/${{SAMPLE_NAME}}.insertion_sites.discarded.bed
        done
        ) 2> {log}
        """