FROM continuumio/miniconda3:latest

# System dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set up conda environment
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && \
    conda clean -a -y
SHELL ["conda", "run", "-n", "samap_env", "/bin/bash", "-c"]

# Set default environment
ENV PATH /opt/conda/envs/samap_env/bin:$PATH

# Optional: install Nextflow (or mount from host)
RUN wget -qO- https://get.nextflow.io | bash && mv nextflow /usr/local/bin/

# Set working directory
WORKDIR /workspace
