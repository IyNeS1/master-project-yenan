FROM ubuntu:22.04

LABEL maintainer="Yenan Wang"
LABEL description="ONT modkit for modification aggregation and validation"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install modkit 0.6.2 from conda-forge/bioconda
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh \
    && /opt/conda/bin/conda install -y -c bioconda -c conda-forge ont-modkit=0.6.2 \
    && /opt/conda/bin/conda clean -afy

ENV PATH="/opt/conda/bin:$PATH"

CMD ["/bin/bash"]
