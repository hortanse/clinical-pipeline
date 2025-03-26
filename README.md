# Clinical Bioinformatics Pipeline

A robust, clinical-grade genomic analysis pipeline for variant detection and interpretation in a healthcare setting.

## Overview

This pipeline implements a comprehensive genomic analysis workflow designed for clinical applications, with a focus on reliability, reproducibility, and compliance with healthcare standards. It processes raw sequencing data through alignment, variant calling, annotation, and clinical interpretation to generate actionable reports for healthcare providers.

## Features

- Sample management and data validation
- Quality control and read preprocessing
- Alignment to reference genome using BWA-MEM
- Variant calling with GATK best practices
- Comprehensive variant annotation using VEP and clinical databases
- Clinical interpretation based on ACMG/AMP guidelines
- Customizable clinical reporting
- EHR integration capabilities
- Audit trail and data provenance tracking
- PHI handling compliance features

## Requirements

- Nextflow (>=22.10.0)
- Java 11 or later
- Docker or Singularity (for containerized execution)
- Reference genome and annotation files
- 16+ CPU cores and 64+ GB RAM recommended for production

### Docker Container Execution

This pipeline is designed to use Docker containers for all third-party tools to ensure reproducibility and ease of deployment. The following containers are used:

- `broadinstitute/gatk:4.3.0.0` - For all GATK-based processes
- `biocontainers/bwa:v0.7.17_cv1` - For read alignment
- `biocontainers/samtools:1.15.1` - For SAM/BAM file processing
- `ensemblorg/ensembl-vep:release_108.1` - For variant annotation
- `nfcore/fastqc:0.11.9` - For quality control
- `staphb/trimgalore:0.6.7` - For read trimming
- `biocontainers/fastp:0.23.2` - For read preprocessing
- `biocontainers/bcftools:1.15.1` - For VCF file manipulation
- `nfcore/snpeff:5.1` - For variant effect prediction

No local installation of these tools is required when using the containerized execution.

### Using Docker with the Pipeline

1. Install Docker:
   ```bash
   # Ubuntu
   curl -fsSL https://get.docker.com -o get-docker.sh
   sudo sh get-docker.sh
   
   # Add your user to the docker group (optional, for running docker without sudo)
   sudo usermod -aG docker $USER
   # Log out and back in for this to take effect
   ```

2. Ensure Docker is running:
   ```bash
   docker --version
   docker run hello-world
   ```

3. Run the pipeline with Docker:
   ```bash
   nextflow run main.nf -profile docker --input samplesheet.csv --output results
   ```

### Container Configuration

Container-specific settings can be adjusted in the `nextflow.config` file:

```nextflow
docker {
    enabled = true
    // Use fixed version tags for reproducibility
    fixRevision = true
    // Remove containers after use
    removeContainers = true
    // Run with elevated privileges (may be needed for some mounts)
    runOptions = '-u $(id -u):$(id -g)'
}
```

### Benefits of Containerized Execution

- **Reproducibility**: Exact same software versions used across all environments
- **Portability**: Run on any system with Docker installed without complex dependencies
- **Isolation**: Tools run in isolated environments to avoid conflicts
- **Version Control**: Specific container versions ensure consistent results
- **Simplified Deployment**: No need to install and configure multiple bioinformatics tools

### Required software (if not using containers)

If you prefer not to use containers, you'll need to install these tools locally:

- BWA
- Samtools
- GATK 4+
- Ensembl VEP
- SnpEff
- BCFtools
- Python 3.8+
- R 4.0+

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/hortanse/clinical-pipeline.git
   cd clinical-pipeline
   ```

2. If not using Docker containers, install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Install Nextflow:
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

4. Run the pipeline:
   ```bash
   ./nextflow run main.nf -profile standard --input samplesheet.csv --output results
   ```

### Docker Installation (Recommended)

For the containerized approach (recommended):

1. Install Docker as described in the Docker section above

2. Build custom Docker images (optional, only needed if modifying the containers):
   ```bash
   # Build main container
   docker build -t clinical-pipeline:latest -f docker/Dockerfile .
   
   # Build VEP container
   docker build -t clinical-pipeline-vep:latest -f docker/vep/Dockerfile docker/vep/
   ```

3. Run the pipeline with Docker profile:
   ```bash
   ./nextflow run main.nf -profile docker --input samplesheet.csv --output results
   ```

This will automatically pull required containers from public repositories and use the local ones you've built.

## Pipeline Structure

```
clinical-pipeline/
├── main.nf                  # Main Nextflow workflow
├── nextflow.config          # Nextflow configuration
├── README.md                # This documentation
├── requirements.txt         # Python dependencies
├── modules/                 # Nextflow process modules
│   ├── ingestion.nf
│   ├── preprocessing.nf
│   ├── alignment.nf
│   ├── variant_calling.nf
│   ├── annotation.nf
│   ├── interpretation.nf
│   ├── reporting.nf
│   └── ehr_integration.nf
├── lib/                     # Python libraries
│   ├── lims/                # LIMS integration
│   ├── variant/             # Variant analysis
│   ├── reporting/           # Report generation
│   └── ehr/                 # EHR integration
├── docker/                  # Docker configurations
│   ├── Dockerfile          # Main Dockerfile
│   ├── base/               # Base image configurations
│   ├── gatk/               # GATK-specific Dockerfile
│   └── vep/                # VEP-specific Dockerfile
├── conf/                    # Configuration profiles
├── assets/                  # Report templates, etc.
└── tests/                   # Test suite
```

## Usage

### Basic Execution

```bash
nextflow run main.nf --input samplesheet.csv --output results
```

### With Docker

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --output results
```

### Advanced Options

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --genome GRCh38 \
  --panel cancer_gene_panel.bed \
  --ehr_integration true \
  --ehr_api_key your_api_key \
  -profile docker,production
```

### Input Format

The input samplesheet should be a CSV file with the following columns:
- sample_id
- patient_id
- fastq_1
- fastq_2
- library
- platform

Example:
```csv
sample_id,patient_id,fastq_1,fastq_2,library,platform
SAMPLE1,PATIENT1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,lib1,ILLUMINA
```

## Configuration Profiles

- `standard`: Local execution (default)
- `docker`: Use Docker containers for all tools
- `singularity`: Use Singularity containers (alternative to Docker)
- `production`: SLURM cluster execution with healthcare compliance settings
- `cloud`: AWS Batch execution
- `high_sensitivity`: Higher sensitivity settings for rare disease applications

Profiles can be combined with comma separation: `-profile docker,production`

## Output Files

The pipeline generates the following output directory structure:

```
results/
├── fastqc/                  # Quality control reports
├── aligned/                 # Aligned BAM files
├── variants/                # Called and filtered variants
│   ├── raw/                 # Raw VCF files
│   ├── filtered/            # Filtered VCF files
│   └── annotated/           # Annotated VCF files
├── interpretation/          # Clinical interpretation files
├── reports/                 # Clinical reports (HTML, PDF, JSON)
│   ├── SAMPLE1.report.html
│   ├── SAMPLE1.report.pdf
│   └── SAMPLE1.report.json
└── pipeline_info/           # Execution reports and logs
```

## Clinical Compliance

This pipeline incorporates features to support clinical genomics operations:

- PHI encryption for sensitive data
- Audit logging for all operations
- Version tracking for all software and databases
- Validation reports
- Data provenance tracking

## Contributors

- Yi-Fan Chou hortanse@gmail.com

## License



## Citation





