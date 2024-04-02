###########################
## FASTQ Trimming and QC ##
###########################

# Trim the R1 and R2 adapters using cutadapt.
# Default TRACE read structure:
#   - R1: [Random Seq (4bp)][SB IR-binding Primer Site][partial SB IR][insertion site][genomic insert]
#   - R2: [Nextera Read 1 adapter][UMI (8 bp)][Tn5 mosaic][genomic insert]
# 1 ) The adapter trimming step checks for:
#   - R1: [Random Seq (4bp)][SB IR-binding Primer Site]
#   - R2: [Nextera Read 1 adapter][UMI (8 bp)][Tn5 mosaic]
# 2) The IR trimming step cehcks for:
#   - R1: [partial SB IR]
# 3) The adapter readthrough step checks for the reverse complement of the R1 and R2 adapters
#   on the 3' end of the R2 and R1 reads respectively.
# 4) The IR readthrough step checks for the reverse complement of the [partial SB IR]
#   on the 3' end of the R2 read.

rule trim_adapters:
    input:
        R1 = fastq_dir / "{sample_name}_R1_001.fastq.gz",
        R2 = fastq_dir / "{sample_name}_R2_001.fastq.gz"
    output:
        # Adapter trimming.
        R1_adapter_trimmed = outdir / 'trimmed/{sample_name}_R1.adapter.trimmed.fastq.gz',
        R2_adapter_trimmed = outdir / 'trimmed/{sample_name}_R2.adapter.trimmed.fastq.gz',
        R1_adapter_untrimmed = outdir / 'trimmed/{sample_name}_R1.adapter.untrimmed.fastq.gz',
        R2_adapter_untrimmed = outdir / 'trimmed/{sample_name}_R2.adapter.untrimmed.fastq.gz',
        adapter_trimming_log = outdir / 'trimmed/{sample_name}.adapter.trimmed.txt',
        # IR trimming.
        R1_ir_trimmed = outdir / 'trimmed/{sample_name}_R1.ir.trimmed.fastq.gz',
        R2_ir_trimmed = outdir / 'trimmed/{sample_name}_R2.ir.trimmed.fastq.gz',
        R1_ir_untrimmed = outdir / 'trimmed/{sample_name}_R1.ir.untrimmed.fastq.gz',
        R2_ir_untrimmed = outdir / 'trimmed/{sample_name}_R2.ir.untrimmed.fastq.gz',
        ir_trimming_log = outdir / 'trimmed/{sample_name}.ir.trimmed.txt',
        # Adapter readthrough trimming.
        R1_adapter_readthrough_trimmed = outdir / 'trimmed/{sample_name}_R1.adapter.readthrough.trimmed.fastq.gz',
        R2_adapter_readthrough_trimmed = outdir / 'trimmed/{sample_name}_R2.adapter.readthrough.trimmed.fastq.gz',
        adapter_readthrough_trimming_log = outdir / 'trimmed/{sample_name}.adapter.readthrough.trimmed.txt',
        # IR readthrough trimming.
        R1_ir_readthrough_trimmed = outdir / 'trimmed/{sample_name}_R1.ir.readthrough.trimmed.fastq.gz',
        R2_ir_readthrough_trimmed = outdir / 'trimmed/{sample_name}_R2.ir.readthrough.trimmed.fastq.gz',
        ir_readthrough_trimming_log = outdir / 'trimmed/{sample_name}.ir.readthrough.trimmed.txt'
    params:
        R1_adapter = R1_adapter,
        R2_adapter = R2_adapter,
        ir_seq = IR_seq,
        trim_adapter_errors = 5,
        trim_ir_errors = 4,
        trim_adapter_readthrough_errors = 4,
        trim_ir_readthrough_errors = 3,
        min_length = 20
    log: outdir / 'logs/{sample_name}.trim_adapters.log'
    benchmark: outdir / 'benchmarks/{sample_name}.trim_adapters.benchmark.txt'
    threads: 6
    conda: '../envs/cutadapt.yaml'
    container: 'docker://quay.io/biocontainers/cutadapt:4.4--py39hf95cd2a_1'
    shell:
        """
        (
        # Adapter trimming.
        cutadapt \
            -j {threads} \
            -g {params.R1_adapter} \
            -G {params.R2_adapter} \
            -e {params.trim_adapter_errors} \
            -m {params.min_length} \
            -o {output.R1_adapter_trimmed} \
            -p {output.R2_adapter_trimmed} \
            --untrimmed-output {output.R1_adapter_untrimmed} \
            --untrimmed-paired-output {output.R2_adapter_untrimmed} \
            {input.R1} \
            {input.R2} \
            > {output.adapter_trimming_log}
        
        # IR trimming.
        cutadapt \
            -j {threads} \
            -g {params.ir_seq} \
            -e {params.trim_ir_errors} \
            -m {params.min_length} \
            -o {output.R1_ir_trimmed} \
            -p {output.R2_ir_trimmed} \
            --untrimmed-output {output.R1_ir_untrimmed} \
            --untrimmed-paired-output {output.R2_ir_untrimmed} \
            {output.R1_adapter_trimmed} \
            {output.R2_adapter_trimmed} \
            > {output.ir_trimming_log}

        # Trim adapter readthrough.
        cutadapt \
            -j {threads} \
            -a $(echo {params.R2_adapter} | sed 's/[^ATGCN]//g' | tr 'ATGCN' 'TACGN' | rev) \
            -A $(echo {params.R1_adapter} | sed 's/[^ATGCN]//g' | tr 'ATGCN' 'TACGN' | rev) \
            -m {params.min_length} \
            -e {params.trim_adapter_readthrough_errors} \
            -o {output.R1_adapter_readthrough_trimmed} \
            -p {output.R2_adapter_readthrough_trimmed} \
            {output.R1_ir_trimmed} \
            {output.R2_ir_trimmed} \
            > {output.adapter_readthrough_trimming_log}
        
        # Trim IR readthrough.
        cutadapt \
            -j {threads} \
            -A $(echo {params.ir_seq} | sed 's/[^ATGCN]//g' | tr 'ATGCN' 'TACGN' | rev) \
            -m {params.min_length} \
            -e {params.trim_ir_readthrough_errors} \
            -o {output.R1_ir_readthrough_trimmed} \
            -p {output.R2_ir_readthrough_trimmed} \
            {output.R1_adapter_readthrough_trimmed} \
            {output.R2_adapter_readthrough_trimmed} \
            > {output.ir_readthrough_trimming_log}        
        ) 2> {log}
        """

# FastQC report of all the reads.
rule fastqc:
    input:
        R1 = fastq_dir / '{sample_name}_R1_001.fastq.gz',
        R2 = fastq_dir / '{sample_name}_R2_001.fastq.gz',
        R1_adapter_trimmed = outdir / 'trimmed/{sample_name}_R1.adapter.trimmed.fastq.gz',
        R2_adapter_trimmed = outdir / 'trimmed/{sample_name}_R2.adapter.trimmed.fastq.gz',
        R1_ir_trimmed = outdir / 'trimmed/{sample_name}_R1.ir.trimmed.fastq.gz',
        R2_ir_trimmed = outdir / 'trimmed/{sample_name}_R2.ir.trimmed.fastq.gz',
        R1_adapter_readthrough_trimmed = outdir / 'trimmed/{sample_name}_R1.adapter.readthrough.trimmed.fastq.gz',
        R2_adapter_readthrough_trimmed = outdir / 'trimmed/{sample_name}_R2.adapter.readthrough.trimmed.fastq.gz',
        R1_ir_readthrough_trimmed = outdir / 'trimmed/{sample_name}_R1.ir.readthrough.trimmed.fastq.gz',
        R2_ir_readthrough_trimmed = outdir / 'trimmed/{sample_name}_R2.ir.readthrough.trimmed.fastq.gz'
    output:
        R1 = outdir / 'fastqc/{sample_name}_R1_001_fastqc.zip',
        R2 = outdir / 'fastqc/{sample_name}_R2_001_fastqc.zip',
        R1_adapter_trimmed = outdir / 'fastqc/{sample_name}_R1.adapter.trimmed_fastqc.zip',
        R2_adapter_trimmed = outdir / 'fastqc/{sample_name}_R2.adapter.trimmed_fastqc.zip',
        R1_ir_trimmed = outdir / 'fastqc/{sample_name}_R1.ir.trimmed_fastqc.zip',
        R2_ir_trimmed = outdir / 'fastqc/{sample_name}_R2.ir.trimmed_fastqc.zip',
        R1_adapter_readthrough_trimmed = outdir / 'fastqc/{sample_name}_R1.adapter.readthrough.trimmed_fastqc.zip',
        R2_adapter_readthrough_trimmed = outdir / 'fastqc/{sample_name}_R2.adapter.readthrough.trimmed_fastqc.zip',
        R1_ir_readthrough_trimmed = outdir / 'fastqc/{sample_name}_R1.ir.readthrough.trimmed_fastqc.zip',
        R2_ir_readthrough_trimmed = outdir / 'fastqc/{sample_name}_R2.ir.readthrough.trimmed_fastqc.zip'
    params:
        outdir = outdir
    log: outdir / 'logs/{sample_name}.fastqc.log'
    benchmark: outdir / 'benchmarks/{sample_name}.fastqc.benchmark.txt'
    threads: 6
    conda: '../envs/fastqc.yaml'
    container: 'docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    shell:
        """
        (
        fastqc \
            -t {threads} \
            -o {params.outdir}/fastqc \
            {input.R1} \
            {input.R2} \
            {input.R1_adapter_trimmed} \
            {input.R2_adapter_trimmed} \
            {input.R1_ir_trimmed} \
            {input.R2_ir_trimmed} \
            {input.R1_adapter_readthrough_trimmed} \
            {input.R2_adapter_readthrough_trimmed} \
            {input.R1_ir_readthrough_trimmed} \
            {input.R2_ir_readthrough_trimmed}
        ) 2> {log}
        """