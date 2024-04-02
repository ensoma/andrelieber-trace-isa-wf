FROM mambaorg/micromamba:1.5.1

RUN \
    micromamba create \
        -y -n is_anno -c conda-forge -c bioconda \
        r-base==4.3.1 \
        r-tidyverse==2.0.0 \
        r-patchwork==1.1.3 \
        r-gridextra==2.3 \
        r-randomcolor==1.1.0.1 \
        r-optparse==1.7.3 \
        r-ggseqlogo==0.1 \
        r-assertthat==0.2.1 \
        bioconductor-biostrings==2.68.1 \
        bioconductor-rtracklayer==1.60.0 \
        bioconductor-genomicfeatures==1.52.1 \
        bioconductor-chipseeker==1.36.0 && \
    micromamba clean -y --all

ENV PATH /opt/conda/envs/is_anno/bin:$PATH