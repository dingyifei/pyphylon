FROM snakemake/snakemake:v8.30.0
WORKDIR /pipeline

# Python scientific stack + BLAST
COPY requirements.txt /tmp/requirements.txt
RUN micromamba install -y -n base -c conda-forge -c bioconda blast less curl wget \
    && micromamba run -n base pip install --no-cache-dir -r /tmp/requirements.txt

# Quarto + TinyTeX for PDF reports
ARG QUARTO_VERSION=1.6.42
RUN curl -fsSL "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    -o /tmp/quarto.deb && dpkg -i /tmp/quarto.deb && rm /tmp/quarto.deb \
    && micromamba run -n base quarto install tinytex --no-prompt

# Entrypoint: install pyphylon from mounted source, then run command
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["snakemake", "--cores", "all", "--software-deployment-method", "apptainer"]
