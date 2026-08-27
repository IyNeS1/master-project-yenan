FROM ubuntu:22.04

LABEL maintainer="Yenan Wang"
LABEL description="ONT modkit for modification aggregation and validation"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    libzstd1 \
    && rm -rf /var/lib/apt/lists/*

# Install modkit 0.6.2 from GitHub releases
RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.6.2/modkit_v0.6.2_u16_x86_64.tar.gz \
    && tar -xzf modkit_v0.6.2_u16_x86_64.tar.gz \
    && mv modkit /usr/local/bin/ \
    && rm modkit_v0.6.2_u16_x86_64.tar.gz

CMD ["/bin/bash"]
