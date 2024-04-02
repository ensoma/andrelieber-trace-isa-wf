FROM mambaorg/micromamba:1.5.1

RUN \
    micromamba create \
        -y -n is_softclip -c conda-forge -c bioconda \
        pysam==0.21.0 \
        samtools==1.17 && \
    micromamba clean -y --all

ENV PATH /opt/conda/envs/is_softclip/bin:$PATH