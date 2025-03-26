FROM ubuntu:22.04

LABEL maintainer="Your Organization <support@example.org>"
LABEL description="Clinical Bioinformatics Pipeline - Base Container"
LABEL version="1.0.0"

# Set noninteractive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install common dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    unzip \
    git \
    ca-certificates \
    build-essential \
    libz-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip \
    default-jre \
    r-base \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up folder structure
RUN mkdir -p /data /references /tools /output

# Create reference symlinks
RUN mkdir -p /references/GRCh38 /references/GRCh37

# Set up Python
RUN pip3 install --no-cache-dir \
    pandas \
    numpy \
    pysam \
    matplotlib \
    biopython \
    pyyaml

# Set working directory
WORKDIR /pipeline

# Set PATH
ENV PATH="/tools:${PATH}"

# Default command
CMD ["bash"] 