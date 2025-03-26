FROM python:3.9-slim

LABEL maintainer="Your Organization <support@example.org>"
LABEL description="Clinical Bioinformatics Pipeline - Main Container"
LABEL version="1.0.0"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    unzip \
    procps \
    git \
    ca-certificates \
    build-essential \
    libz-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /pipeline

# Install Python dependencies for the pipeline
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Add Python library path to PYTHONPATH
ENV PYTHONPATH="/pipeline:${PYTHONPATH}"

# Default command
CMD ["bash"] 