FROM nfcore/base:1.13.2
LABEL authors="Leon Bichmann" \
      description="Docker image containing all software requirements for the nf-core/diaproteomics pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-diaproteomics-dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-diaproteomics-dev > nf-core-diaproteomics-dev.yml
