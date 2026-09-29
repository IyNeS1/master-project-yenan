FROM ubuntu:22.04

ARG DORADO_VERSION=1.3.1

LABEL org.opencontainers.image.source="https://github.com/IyNeS1/master-project-yenan"
LABEL org.opencontainers.image.description="Dorado ${DORADO_VERSION} basecaller for ONT direct RNA (RNA004)"

# procps is required by Nextflow for task metrics (ps command)
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates curl procps && \
    rm -rf /var/lib/apt/lists/*

# Download the official Dorado release; it bundles its own CUDA/libtorch libraries
RUN curl -fsSL "https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz" \
    | tar -xz -C /opt

# Use PATH instead of a symlink so Dorado resolves its bundled lib/ directory correctly
ENV PATH="/opt/dorado-${DORADO_VERSION}-linux-x64/bin:${PATH}"
