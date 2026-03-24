FROM condaforge/miniforge3:latest
WORKDIR /pipeline

# Layer 1 — system packages
RUN DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/share/zoneinfo/UTC /etc/localtime

# Layer 2 — conda base env (snakemake + blast + ML packages)
RUN rm -rf /opt/conda/pkgs/cache/*.json \
    && mamba install -y -n base -c conda-forge -c bioconda \
      snakemake blast=2.16.0 less curl wget \
      hdbscan cxx-compiler \
    && pip install --no-cache-dir fastcluster \
    && mamba clean -afy

# Layer 3 — pinned Python science/viz stack
COPY requirements-core.txt /tmp/requirements-core.txt
RUN pip install --no-cache-dir -r /tmp/requirements-core.txt

# Layer 4 — marimo, LSP, dev tools
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

EXPOSE 2718

COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["snakemake", "--cores", "all", "--sdm", "conda"]
