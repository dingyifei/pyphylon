FROM condaforge/mambaforge:latest
WORKDIR /pipeline

# Python scientific stack + bioinformatics tools + TeX for PDF reports
COPY requirements.txt /tmp/requirements.txt
RUN mamba install -y -n base -c conda-forge -c bioconda \
      snakemake apptainer blast less curl wget \
      hdbscan cxx-compiler texlive-core \
    && pip install --no-cache-dir fastcluster -r /tmp/requirements.txt \
    && mamba clean -afy \
    && DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/share/zoneinfo/UTC /etc/localtime

# Quarto (arch-aware download)
ARG QUARTO_VERSION=1.6.42
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "aarch64" ] || [ "$ARCH" = "arm64" ]; then \
      QUARTO_ARCH="linux-arm64"; \
    else \
      QUARTO_ARCH="linux-amd64"; \
    fi && \
    curl -fsSL "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-${QUARTO_ARCH}.deb" \
      -o /tmp/quarto.deb && dpkg -i /tmp/quarto.deb && rm /tmp/quarto.deb

# Entrypoint: install pyphylon from mounted source, then run command
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["snakemake", "--cores", "all", "--software-deployment-method", "apptainer"]
