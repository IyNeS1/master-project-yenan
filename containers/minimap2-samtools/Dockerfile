FROM ubuntu:22.04

LABEL maintainer="Yenan Wang"
LABEL description="minimap2 and samtools for ONT direct RNA sequencing alignment"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    gcc \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libncurses5-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install minimap2 2.30
RUN wget https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30_x64-linux.tar.bz2 \
    && tar -xjf minimap2-2.30_x64-linux.tar.bz2 \
    && cp minimap2-2.30_x64-linux/minimap2 /usr/local/bin/ \
    && rm -rf minimap2-2.30*

# Install samtools 1.23.1
RUN wget https://github.com/samtools/samtools/releases/download/1.23.1/samtools-1.23.1.tar.bz2 \
    && tar -xjf samtools-1.23.1.tar.bz2 \
    && cd samtools-1.23.1 \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf samtools-1.23.1*

CMD ["/bin/bash"]
