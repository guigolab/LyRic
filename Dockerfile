FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="ea274b87e1e23c8012338f83b0549a1d778db448db86f5239f446032503c09c0"
LABEL org.opencontainers.image.source="https://github.com/guigolab/LyRic"

# Conda environment:
#   name: R_env
RUN mkdir -p /conda-envs/88953de6d96581fc47f6cfb8b1181c40
COPY envs/R_env.yml /conda-envs/88953de6d96581fc47f6cfb8b1181c40/environment.yaml

# Conda environment:
#   name: Redtools_env
RUN mkdir -p /conda-envs/455e8d2b42ec6004358e1233391fd7b8
COPY envs/Redtools_env.yml /conda-envs/455e8d2b42ec6004358e1233391fd7b8/environment.yaml

# Conda environment:
#   name: gffcompare_env
RUN mkdir -p /conda-envs/4754a1a6e178d2e5c343134e5a2ccd57
COPY envs/gffcompare_env.yml /conda-envs/4754a1a6e178d2e5c343134e5a2ccd57/environment.yaml

# Conda environment:
#   name: minimap2_env
RUN mkdir -p /conda-envs/cee76cc0c4070e4b8a385824f7151b2c
COPY envs/minimap2_env.yml /conda-envs/cee76cc0c4070e4b8a385824f7151b2c/environment.yaml

# Conda environment:
#   name: perl_env
RUN mkdir -p /conda-envs/d269bede2efc1304ae0d31b80dd637cf
COPY envs/perl_env.yml /conda-envs/d269bede2efc1304ae0d31b80dd637cf/environment.yaml

# Conda environment:
#   name: qualimap_env
RUN mkdir -p /conda-envs/6c2a3938a71045f4f82e72c2749a9a11
COPY envs/qualimap_env.yml /conda-envs/6c2a3938a71045f4f82e72c2749a9a11/environment.yaml

# Conda environment:
#   name: ucsc_env
RUN mkdir -p /conda-envs/26b70ccebb2c14273b0a42abc1e376b4
COPY envs/ucsc_env.yml /conda-envs/26b70ccebb2c14273b0a42abc1e376b4/environment.yaml

# Conda environment:
#   name: xtools_env
RUN mkdir -p /conda-envs/fc27316f1f50daf859adba261be14d3c
COPY envs/xtools_env.yml /conda-envs/fc27316f1f50daf859adba261be14d3c/environment.yaml

# Generate conda environments
RUN mamba env create --prefix /conda-envs/88953de6d96581fc47f6cfb8b1181c40 --file /conda-envs/88953de6d96581fc47f6cfb8b1181c40/environment.yaml && \
    mamba env create --prefix /conda-envs/455e8d2b42ec6004358e1233391fd7b8 --file /conda-envs/455e8d2b42ec6004358e1233391fd7b8/environment.yaml && \
    mamba env create --prefix /conda-envs/4754a1a6e178d2e5c343134e5a2ccd57 --file /conda-envs/4754a1a6e178d2e5c343134e5a2ccd57/environment.yaml && \
    mamba env create --prefix /conda-envs/cee76cc0c4070e4b8a385824f7151b2c --file /conda-envs/cee76cc0c4070e4b8a385824f7151b2c/environment.yaml && \
    mamba env create --prefix /conda-envs/d269bede2efc1304ae0d31b80dd637cf --file /conda-envs/d269bede2efc1304ae0d31b80dd637cf/environment.yaml && \
    mamba env create --prefix /conda-envs/6c2a3938a71045f4f82e72c2749a9a11 --file /conda-envs/6c2a3938a71045f4f82e72c2749a9a11/environment.yaml && \
    mamba env create --prefix /conda-envs/26b70ccebb2c14273b0a42abc1e376b4 --file /conda-envs/26b70ccebb2c14273b0a42abc1e376b4/environment.yaml && \
    mamba env create --prefix /conda-envs/fc27316f1f50daf859adba261be14d3c --file /conda-envs/fc27316f1f50daf859adba261be14d3c/environment.yaml && \
    mamba clean --all -y

RUN mamba install conda-forge::gawk \
    && mamba clean --all -y

RUN apt update && \
    apt install uuid-runtime && \
    apt clean
