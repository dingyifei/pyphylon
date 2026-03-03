FROM condaforge/mambaforge:latest
WORKDIR /pipeline

COPY requirements.txt /tmp/requirements.txt
RUN mamba install -y -n base -c conda-forge -c bioconda \
      snakemake apptainer blast less curl wget \
      hdbscan cxx-compiler \
    && pip install --no-cache-dir fastcluster -r /tmp/requirements.txt \
    && mamba clean -afy \
    && DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/share/zoneinfo/UTC /etc/localtime

EXPOSE 2718

COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["snakemake", "--cores", "all", "--software-deployment-method", "apptainer"]
