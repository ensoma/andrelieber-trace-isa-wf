#############
## MultiQC ##
#############

# Final MultiQC report.
rule final_multiqc_report:
    input:
        # FastQC
        fastqc_R1_raw = [outdir / f"fastqc/{sample}_R1_001_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R2_raw = [outdir / f"fastqc/{sample}_R2_001_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R1_adapter_trimmed = [outdir / f"fastqc/{sample}_R1.adapter.trimmed_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R2_adapter_trimmed = [outdir / f"fastqc/{sample}_R2.adapter.trimmed_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R1_ir_trimmed = [outdir / f"fastqc/{sample}_R1.ir.trimmed_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R2_ir_trimmed = [outdir / f"fastqc/{sample}_R2.ir.trimmed_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R1_adapter_readthrough_trimmed = [outdir / f"fastqc/{sample}_R1.adapter.readthrough.trimmed_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R2_adapter_readthrough_trimmed = [outdir / f"fastqc/{sample}_R2.adapter.readthrough.trimmed_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R1_ir_readthrough_trimmed = [outdir / f"fastqc/{sample}_R1.ir.readthrough.trimmed_fastqc.zip" for sample in sample_map.keys()],
        fastqc_R2_ir_readthrough_trimmed = [outdir / f"fastqc/{sample}_R2.ir.readthrough.trimmed_fastqc.zip" for sample in sample_map.keys()],
        # Cutadapt trimmed reads.
        adapter_trimmed_log = [outdir / f"trimmed/{sample}.adapter.trimmed.txt" for sample in sample_map.keys()],
        ir_trimmed_log = [outdir / f"trimmed/{sample}.ir.trimmed.txt" for sample in sample_map.keys()],
        adapter_readthrough_trimmed_log = [outdir / f"trimmed/{sample}.adapter.readthrough.trimmed.txt" for sample in sample_map.keys()],
        ir_readthrough_trimmed_log = [outdir / f"trimmed/{sample}.ir.readthrough.trimmed.txt" for sample in sample_map.keys()],
        # Alignment stats.
        sam_flagstat = [outdir / f"alignments/{sample}.aligned.sam.flagstat" for sample in sample_map.keys()],
        sam_idxstat = [outdir / f"alignments/{sample}.aligned.sam.idxstat" for sample in sample_map.keys()],
        filtered_bam_flagstat = [outdir / f"alignments/{sample}.filtered.bam.flagstat" for sample in sample_map.keys()],
        filtered_bam_idxstat = [outdir / f"alignments/{sample}.filtered.bam.idxstat" for sample in sample_map.keys()],
        deduped_bam_flagstat = [outdir / f"alignments/{sample}.deduped.bam.flagstat" for sample in sample_map.keys()],
        deduped_bam_idxstat = [outdir / f"alignments/{sample}.deduped.bam.idxstat" for sample in sample_map.keys()],
        softclip_filtered_bam_flagstat = [outdir / f"alignments/{sample}.softclip_removed.bam.flagstat" for sample in sample_map.keys()],
        softclip_filtered_bam_idxstat = [outdir / f"alignments/{sample}.softclip_removed.bam.idxstat" for sample in sample_map.keys()],
        # Aggregate stats.
        aggregate_stats = outdir / 'stats/aggregate_stats.tsv',
        aggregate_mapping_alignment_stats = outdir / 'stats/aggregate.trimming_and_mapping_stats.tsv',
        # Seqlogo.
        seqlogo = outdir / 'seqlogo/is_seqlogo.png',
        # IS fractions stacked barplot.
        is_fractions = outdir / 'insertion_sites/is_fractions_plot.png',
        # Insertion site annotation plots.
        feature_dist_table = outdir / 'annotated_insertion_sites/is_feature_distribution_table.tsv',
        promoter_dist_plot = outdir / 'annotated_insertion_sites/promoter_distribution_plot.png',
        top_is_table = outdir / 'annotated_insertion_sites/top_insertion_sites_table.tsv',
        # Blacklist aggregated stats.
        blacklisted = outdir / 'insertion_sites/discarded_is_report.tsv'
    output:
        outdir / 'multiqc/multiqc_report.html'
    params:
        outdir = outdir / 'multiqc',
        multiqc_config = 'config/multiqc_config.yaml'
    log: outdir / 'logs/final_multiqc_report.log'
    benchmark: outdir / 'benchmarks/final_multiqc_report.benchmark.txt'
    conda: '../envs/multiqc.yaml'
    container: 'docker://quay.io/biocontainers/multiqc:1.16--pyhdfd78af_0'
    shell:
        """
        (
        multiqc \
            -o {params.outdir} \
            -f \
            -c {params.multiqc_config} \
            {input.fastqc_R1_raw} \
            {input.fastqc_R2_raw} \
            {input.fastqc_R1_adapter_trimmed} \
            {input.fastqc_R2_adapter_trimmed} \
            {input.fastqc_R1_ir_trimmed} \
            {input.fastqc_R2_ir_trimmed} \
            {input.fastqc_R1_adapter_readthrough_trimmed} \
            {input.fastqc_R2_adapter_readthrough_trimmed} \
            {input.fastqc_R1_ir_readthrough_trimmed} \
            {input.fastqc_R2_ir_readthrough_trimmed} \
            {input.adapter_trimmed_log} \
            {input.ir_trimmed_log} \
            {input.adapter_readthrough_trimmed_log} \
            {input.ir_readthrough_trimmed_log} \
            {input.sam_flagstat} \
            {input.sam_idxstat} \
            {input.filtered_bam_flagstat} \
            {input.filtered_bam_idxstat} \
            {input.deduped_bam_flagstat} \
            {input.deduped_bam_idxstat} \
            {input.softclip_filtered_bam_flagstat} \
            {input.softclip_filtered_bam_idxstat} \
            {input.aggregate_stats} \
            {input.aggregate_mapping_alignment_stats} \
            {input.seqlogo} \
            {input.is_fractions} \
            {input.feature_dist_table} \
            {input.promoter_dist_plot} \
            {input.top_is_table} \
            {input.blacklisted}
        ) 2> {log}
        """