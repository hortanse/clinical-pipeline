FROM broadinstitute/gatk:4.3.0.0

LABEL maintainer="Your Organization <support@example.org>"
LABEL description="Clinical Bioinformatics Pipeline - GATK Container"
LABEL version="1.0.0"

# Add custom scripts for the pipeline
COPY scripts/ /usr/local/bin/

# Make scripts executable
RUN chmod +x /usr/local/bin/*.sh

# Create reference directory and set up structure
RUN mkdir -p /references/GRCh38 /references/GRCh37

# Set working directory
WORKDIR /data

# Default command
CMD ["gatk", "--help"] 