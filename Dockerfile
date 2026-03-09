FROM condaforge/mambaforge:latest
WORKDIR /pipeline

# Layer 1 — system packages (almost never changes)
RUN DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/share/zoneinfo/UTC /etc/localtime

# Layer 2 — conda bioinformatics tools (very stable)
RUN mamba install -y -n base -c conda-forge -c bioconda \
      snakemake apptainer blast cd-hit less curl wget \
    && mamba clean -afy

# Layer 3 — conda compiled ML packages (stable)
RUN mamba install -y -n base -c conda-forge \
      hdbscan cxx-compiler \
    && pip install --no-cache-dir fastcluster \
    && mamba clean -afy

# Layer 4 — pinned Python science/viz stack (changes infrequently)
COPY requirements-core.txt /tmp/requirements-core.txt
RUN pip install --no-cache-dir -r /tmp/requirements-core.txt

# Layer 5 — marimo, LSP, dev tools (changes most often)
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

EXPOSE 2718

COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["snakemake", "--cores", "all", "--software-deployment-method", "apptainer"]
