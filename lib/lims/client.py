#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
LIMS API Client for Clinical Bioinformatics Pipeline

This module provides a client for interacting with the Laboratory Information 
Management System (LIMS) to retrieve sample information and metadata.
"""

import os
import sys
import json
import time
import logging
import requests
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("LIMS_CLIENT")

class LIMSClient:
    """Client for interacting with the LIMS API."""
    
    def __init__(self, api_url: str, api_key: str, verify_ssl: bool = True):
        """
        Initialize the LIMS client.
        
        Args:
            api_url: Base URL for the LIMS API.
            api_key: API key for authentication.
            verify_ssl: Whether to verify SSL certificates.
        """
        self.api_url = api_url
        self.api_key = api_key
        self.verify_ssl = verify_ssl
        self.session = self._create_session()
        
    def _create_session(self) -> requests.Session:
        """Create and configure a requests session."""
        session = requests.Session()
        session.headers.update({
            'Authorization': f'Bearer {self.api_key}',
            'Content-Type': 'application/json',
            'Accept': 'application/json'
        })
        return session
    
    def get_sample_metadata(self, sample_id: str) -> Dict:
        """
        Retrieve metadata for a specific sample.
        
        Args:
            sample_id: The unique identifier for the sample.
            
        Returns:
            Dict containing sample metadata.
            
        Raises:
            requests.HTTPError: If the API request fails.
        """
        endpoint = f"{self.api_url}/samples/{sample_id}"
        response = self.session.get(endpoint, verify=self.verify_ssl)
        response.raise_for_status()
        return response.json()
    
    def get_samples_from_batch(self, batch_id: str) -> List[Dict]:
        """
        Retrieve all samples associated with a batch.
        
        Args:
            batch_id: The batch identifier.
            
        Returns:
            List of sample metadata dictionaries.
            
        Raises:
            requests.HTTPError: If the API request fails.
        """
        endpoint = f"{self.api_url}/batches/{batch_id}/samples"
        response = self.session.get(endpoint, verify=self.verify_ssl)
        response.raise_for_status()
        return response.json()
    
    def download_fastq(self, sample_id: str, output_dir: Union[str, Path]) -> Tuple[Path, Path]:
        """
        Download FASTQ files for a sample.
        
        Args:
            sample_id: The unique identifier for the sample.
            output_dir: Directory where FASTQ files will be saved.
            
        Returns:
            Tuple containing paths to the downloaded R1 and R2 FASTQ files.
            
        Raises:
            requests.HTTPError: If the API request fails.
            IOError: If file download or writing fails.
        """
        # Ensure output directory exists
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Get download URLs
        endpoint = f"{self.api_url}/samples/{sample_id}/fastq"
        response = self.session.get(endpoint, verify=self.verify_ssl)
        response.raise_for_status()
        fastq_data = response.json()
        
        # Download R1 file
        r1_url = fastq_data.get('r1_url')
        r1_path = output_dir / f"{sample_id}_R1.fastq.gz"
        self._download_file(r1_url, r1_path)
        logger.info(f"Downloaded R1 FASTQ to {r1_path}")
        
        # Download R2 file
        r2_url = fastq_data.get('r2_url')
        r2_path = output_dir / f"{sample_id}_R2.fastq.gz"
        self._download_file(r2_url, r2_path)
        logger.info(f"Downloaded R2 FASTQ to {r2_path}")
        
        return r1_path, r2_path
    
    def _download_file(self, url: str, path: Path) -> None:
        """
        Download a file from a URL to a specified path.
        
        Args:
            url: URL to download from.
            path: Path where file will be saved.
            
        Raises:
            requests.HTTPError: If the download request fails.
            IOError: If writing to file fails.
        """
        with self.session.get(url, stream=True, verify=self.verify_ssl) as resp:
            resp.raise_for_status()
            with open(path, 'wb') as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
    
    def update_sample_status(self, sample_id: str, status: str, message: Optional[str] = None) -> Dict:
        """
        Update the status of a sample in the LIMS system.
        
        Args:
            sample_id: The unique identifier for the sample.
            status: New status for the sample (e.g., 'PROCESSING', 'COMPLETED', 'FAILED').
            message: Optional message providing more details about the status.
            
        Returns:
            Dict containing the API response.
            
        Raises:
            requests.HTTPError: If the API request fails.
        """
        endpoint = f"{self.api_url}/samples/{sample_id}/status"
        data = {
            'status': status,
            'timestamp': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        }
        
        if message:
            data['message'] = message
        
        response = self.session.post(endpoint, json=data, verify=self.verify_ssl)
        response.raise_for_status()
        return response.json()


def parse_samplesheet(samplesheet_path: Union[str, Path]) -> List[Dict]:
    """
    Parse a CSV samplesheet file used to specify samples for processing.
    
    Expected CSV format:
    sample_id,patient_id,batch_id,panel,priority,additional_metadata
    
    Args:
        samplesheet_path: Path to the CSV samplesheet.
        
    Returns:
        List of dictionaries containing sample information.
        
    Raises:
        FileNotFoundError: If samplesheet file does not exist.
        ValueError: If samplesheet format is invalid.
    """
    import csv
    from pathlib import Path
    
    samplesheet_path = Path(samplesheet_path)
    if not samplesheet_path.exists():
        raise FileNotFoundError(f"Samplesheet file not found: {samplesheet_path}")
    
    samples = []
    required_headers = {'sample_id', 'patient_id', 'batch_id'}
    
    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)
        headers = set(reader.fieldnames)
        
        if not required_headers.issubset(headers):
            missing = required_headers - headers
            raise ValueError(f"Samplesheet missing required headers: {', '.join(missing)}")
        
        for row in reader:
            samples.append(row)
    
    return samples


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="LIMS interaction utility")
    parser.add_argument('--samplesheet', required=True, help='Path to sample sheet CSV')
    parser.add_argument('--lims_api_key', required=True, help='API key for LIMS authentication')
    parser.add_argument('--lims_api_url', default=os.environ.get('LIMS_API_URL', 'https://lims.example.org/api/v1'), 
                        help='LIMS API URL')
    parser.add_argument('--output', required=True, help='Output JSON file for sample metadata')
    parser.add_argument('--fastq_dir', required=True, help='Directory to download FASTQ files')
    parser.add_argument('--verify_ssl', action='store_true', default=True, help='Verify SSL certificates')
    
    args = parser.parse_args()
    
    try:
        # Initialize LIMS client
        client = LIMSClient(args.lims_api_url, args.lims_api_key, args.verify_ssl)
        
        # Parse samplesheet
        samples = parse_samplesheet(args.samplesheet)
        
        # Process each sample
        processed_samples = []
        for sample_info in samples:
            sample_id = sample_info['sample_id']
            logger.info(f"Processing sample: {sample_id}")
            
            try:
                # Get detailed metadata
                metadata = client.get_sample_metadata(sample_id)
                
                # Update sample status to PROCESSING
                client.update_sample_status(sample_id, "PROCESSING", "Pipeline processing started")
                
                # Download FASTQ files
                r1_path, r2_path = client.download_fastq(sample_id, args.fastq_dir)
                
                # Add file paths to metadata
                metadata['fastq_r1'] = str(r1_path)
                metadata['fastq_r2'] = str(r2_path)
                
                processed_samples.append(metadata)
                
            except Exception as e:
                logger.error(f"Error processing sample {sample_id}: {str(e)}")
                client.update_sample_status(sample_id, "FAILED", f"Pipeline initialization failed: {str(e)}")
                raise
        
        # Write compiled sample data to output JSON
        with open(args.output, 'w') as f:
            json.dump(processed_samples, f, indent=2)
        
        logger.info(f"Successfully processed {len(processed_samples)} samples")
        
    except Exception as e:
        logger.error(f"Pipeline initialization failed: {str(e)}")
        sys.exit(1) 