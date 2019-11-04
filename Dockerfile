FROM nfcore/base:1.7
LABEL authors="Leon Bichmann" \
      description="Docker image containing all requirements for nf-core/diaproteomics pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-diaproteomics-1.0dev/bin:$PATH
