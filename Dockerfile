# unified toolbox for your pipeline
FROM mambaorg/micromamba:1.5.8

SHELL ["/bin/bash", "-lc"]
USER root

# Base OS tools (note: use gawk, not 'awk')
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      bash coreutils findutils grep sed gawk perl \
      curl wget unzip pigz procps ca-certificates tzdata locales tini \
      tar gzip xz-utils \
  && rm -rf /var/lib/apt/lists/*

# (optional) quiet locale
RUN locale-gen en_US.UTF-8 || true
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# Conda env with all bio tools
COPY env/base.yml /tmp/base.yml
RUN micromamba install -y -n base -f /tmp/base.yml && micromamba clean -ya

# Make conda tools available
ENV PATH=/opt/conda/bin:$PATH

# Writable dirs for arbitrary UID (Nextflow -u)
RUN mkdir -p /data /ref /opt/bin && chmod -R a+rwx /data /ref /tmp

ENTRYPOINT ["/usr/bin/tini","-g","--"]
CMD ["/bin/bash"]

