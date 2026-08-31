FROM python:3.13-slim

LABEL maintainer="Yenan Wang"
LABEL description="Python environment for calibration sweep and summary scripts"

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhts-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir pandas==3.0.2 matplotlib pysam scipy

CMD ["/bin/bash"]
