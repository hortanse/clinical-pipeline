FROM ensemblorg/ensembl-vep:release_108.1

LABEL maintainer="Your Organization <support@example.org>"
LABEL description="Clinical Bioinformatics Pipeline - VEP Container"
LABEL version="1.0.0"

# Install additional plugins and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libdbi-perl \
    libarchive-zip-perl \
    curl \
    unzip \
    tabix \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install additional VEP plugins
RUN mkdir -p /plugins && \
    cd /plugins && \
    curl -L -O https://github.com/Ensembl/VEP_plugins/archive/release/108.zip && \
    unzip 108.zip && \
    rm 108.zip

# Copy plugin configuration
COPY vep_plugins.conf /opt/vep/.vep/plugins.conf

# Set up reference structure
RUN mkdir -p /references/GRCh38 /references/GRCh37

# Set working directory
WORKDIR /data

# Default command
CMD ["vep", "--help"] 