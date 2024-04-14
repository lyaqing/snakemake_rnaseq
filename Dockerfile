FROM continuumio/miniconda:4.7.12

WORKDIR /home/snakemake/

COPY ["./workflow/envs/environment.yaml", "./"]

# mamba is a faster C++ re-implemenbtation of conda
# name of the environment is rnaseq
RUN conda install -c conda-forge mamba --yes \
  && mamba env create -f ./workflow/envs/environment.yaml \
  && conda clean --all

RUN echo "source activate rnaseq" > ~/.bashrc
ENV PATH /opt/conda/envs/rnaseq/bin:$PATH

ENTRYPOINT ["snakemake"]
