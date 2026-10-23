# LEGACY/UNVERIFIED: this image uses environment.yml/micromamba, not pixi.lock.
# Do not treat it as a reproducible deployment; use the verified Pixi host path.
# PfVarX container image (RESEARCH USE ONLY).
#
# Build:
#   docker build -t pfvarx:0.9.0 -f container/Dockerfile .
# Run:
#   docker run --rm -v "$PWD:/data" pfvarx:0.9.0 status
#   docker run --rm -v "$PWD:/data" pfvarx:0.9.0 demo --config /opt/pfvarx/config/profiles/wgs.yaml --out /data/demo

FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="PfVarX"
LABEL org.opencontainers.image.version="0.9.0"
LABEL org.opencontainers.image.description="P. falciparum SNP/INDEL analysis platform (RESEARCH USE ONLY)"

USER root
WORKDIR /opt/pfvarx

# Install the conda environment first (cached layer).
COPY environment.yml /opt/pfvarx/environment.yml
RUN micromamba install -y -n base -f /opt/pfvarx/environment.yml && \
    micromamba clean --all --yes

# Copy the source and install the package.
COPY . /opt/pfvarx
RUN /opt/conda/bin/pip install --no-deps -e /opt/pfvarx && \
    /opt/conda/bin/python -c "import pfvarx; print(pfvarx.status_banner())"

ENV PATH=/opt/conda/bin:$PATH
ENV PFVARX_HOME=/opt/pfvarx

WORKDIR /data
ENTRYPOINT ["/opt/conda/bin/pfvarx"]
CMD ["status"]
