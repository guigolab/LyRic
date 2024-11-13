FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="b9671748f7ff0f3a437a4c751b3a04da8dc6411f6291925400b41d49dfdde77c"
LABEL org.opencontainers.image.source="https://github.com/guigolab/LyRic"

# Conda environment:
#   name: minimap2_env
RUN mkdir -p /conda-envs/cee76cc0c4070e4b8a385824f7151b2c
COPY workflow/envs/minimap2_env.yml /conda-envs/cee76cc0c4070e4b8a385824f7151b2c/environment.yaml

# Conda environment:
#   name: perl_env
RUN mkdir -p /conda-envs/d269bede2efc1304ae0d31b80dd637cf
COPY workflow/envs/perl_env.yml /conda-envs/d269bede2efc1304ae0d31b80dd637cf/environment.yaml

# Conda environment:
#   name: ucsc_env
RUN mkdir -p /conda-envs/26b70ccebb2c14273b0a42abc1e376b4
COPY workflow/envs/ucsc_env.yml /conda-envs/26b70ccebb2c14273b0a42abc1e376b4/environment.yaml

# Conda environment:
#   name: xtools_env
RUN mkdir -p /conda-envs/fc27316f1f50daf859adba261be14d3c
COPY workflow/envs/xtools_env.yml /conda-envs/fc27316f1f50daf859adba261be14d3c/environment.yaml

# Generate conda environments
RUN mamba env create --prefix /conda-envs/cee76cc0c4070e4b8a385824f7151b2c --file /conda-envs/cee76cc0c4070e4b8a385824f7151b2c/environment.yaml && \
    mamba env create --prefix /conda-envs/d269bede2efc1304ae0d31b80dd637cf --file /conda-envs/d269bede2efc1304ae0d31b80dd637cf/environment.yaml && \
    mamba env create --prefix /conda-envs/26b70ccebb2c14273b0a42abc1e376b4 --file /conda-envs/26b70ccebb2c14273b0a42abc1e376b4/environment.yaml && \
    mamba env create --prefix /conda-envs/fc27316f1f50daf859adba261be14d3c --file /conda-envs/fc27316f1f50daf859adba261be14d3c/environment.yaml && \
    mamba clean --all -y

RUN mamba install conda-forge::gawk \
    && mamba clean --all -y

RUN apt update && \
    apt install uuid-runtime && \
    apt clean
