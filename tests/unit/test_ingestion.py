#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the LIMS API client module.
"""

import os
import json
import tempfile
import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path

import requests
import pytest

# Import the module under test
from clinical_pipeline.lib.lims.client import LIMSClient, parse_samplesheet

# Fixtures for testing
SAMPLE_METADATA = {
    "sample_id": "SAMPLE001",
    "patient_id": "PATIENT001",
    "batch_id": "BATCH001",
    "specimen_type": "blood",
    "collection_date": "2023-01-15",
    "received_date": "2023-01-16",
    "panel": "cancer_panel_v1",
    "status": "RECEIVED"
}

SAMPLE_BATCH = [
    SAMPLE_METADATA,
    {
        "sample_id": "SAMPLE002",
        "patient_id": "PATIENT002",
        "batch_id": "BATCH001",
        "specimen_type": "tissue",
        "collection_date": "2023-01-15",
        "received_date": "2023-01-16",
        "panel": "cancer_panel_v1",
        "status": "RECEIVED"
    }
]

FASTQ_URLS = {
    "r1_url": "https://lims.example.org/api/v1/files/SAMPLE001_R1.fastq.gz",
    "r2_url": "https://lims.example.org/api/v1/files/SAMPLE001_R2.fastq.gz"
}


class TestLIMSClient(unittest.TestCase):
    """Test cases for the LIMS API client."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.api_url = "https://lims.example.org/api/v1"
        self.api_key = "test_api_key"
        self.client = LIMSClient(self.api_url, self.api_key, verify_ssl=False)
    
    @patch('requests.Session.get')
    def test_get_sample_metadata(self, mock_get):
        """Test retrieving sample metadata."""
        # Setup mock response
        mock_response = MagicMock()
        mock_response.json.return_value = SAMPLE_METADATA
        mock_get.return_value = mock_response
        
        # Call the method under test
        result = self.client.get_sample_metadata("SAMPLE001")
        
        # Verify results
        mock_get.assert_called_once_with(
            "https://lims.example.org/api/v1/samples/SAMPLE001",
            verify=False
        )
        self.assertEqual(result, SAMPLE_METADATA)
    
    @patch('requests.Session.get')
    def test_get_samples_from_batch(self, mock_get):
        """Test retrieving samples from a batch."""
        # Setup mock response
        mock_response = MagicMock()
        mock_response.json.return_value = SAMPLE_BATCH
        mock_get.return_value = mock_response
        
        # Call the method under test
        result = self.client.get_samples_from_batch("BATCH001")
        
        # Verify results
        mock_get.assert_called_once_with(
            "https://lims.example.org/api/v1/batches/BATCH001/samples",
            verify=False
        )
        self.assertEqual(result, SAMPLE_BATCH)
    
    @patch('requests.Session.get')
    @patch('requests.Session.post')
    def test_update_sample_status(self, mock_post, mock_get):
        """Test updating sample status."""
        # Setup mock response
        mock_response = MagicMock()
        mock_response.json.return_value = {"status": "updated"}
        mock_post.return_value = mock_response
        
        # Call the method under test
        result = self.client.update_sample_status(
            "SAMPLE001", "PROCESSING", "Sample processing started"
        )
        
        # Verify results
        self.assertIn("status", result)
        self.assertEqual(result["status"], "updated")
        
        # Verify the correct endpoint was called
        mock_post.assert_called_once()
        call_args = mock_post.call_args
        self.assertEqual(
            call_args[0][0],
            "https://lims.example.org/api/v1/samples/SAMPLE001/status"
        )
        
        # Verify the correct JSON payload was sent
        json_data = call_args[1]["json"]
        self.assertEqual(json_data["status"], "PROCESSING")
        self.assertEqual(json_data["message"], "Sample processing started")
    
    @patch('requests.Session.get')
    def test_download_fastq(self, mock_get):
        """Test downloading FASTQ files."""
        # Setup mock responses
        fastq_info_response = MagicMock()
        fastq_info_response.json.return_value = FASTQ_URLS
        
        file_response = MagicMock()
        file_response.iter_content.return_value = [b"FASTQ content"]
        
        # Configure the mock to return different responses
        mock_get.side_effect = [
            fastq_info_response,  # First call gets the URLs
            file_response,        # Second call downloads R1
            file_response         # Third call downloads R2
        ]
        
        # Create a temporary directory for downloads
        with tempfile.TemporaryDirectory() as tmpdir:
            # Call the method under test
            r1_path, r2_path = self.client.download_fastq("SAMPLE001", tmpdir)
            
            # Verify the files were downloaded
            self.assertTrue(os.path.exists(r1_path))
            self.assertTrue(os.path.exists(r2_path))
            
            # Verify the correct endpoints were called
            self.assertEqual(mock_get.call_count, 3)
            calls = mock_get.call_args_list
            self.assertEqual(
                calls[0][0][0],
                "https://lims.example.org/api/v1/samples/SAMPLE001/fastq"
            )


class TestSamplesheetParsing:
    """Test cases for samplesheet parsing using pytest style."""
    
    def test_parse_valid_samplesheet(self, tmp_path):
        """Test parsing a valid samplesheet."""
        # Create a sample CSV file
        csv_content = "sample_id,patient_id,batch_id,panel\n" \
                      "SAMPLE001,PATIENT001,BATCH001,cancer_panel_v1\n" \
                      "SAMPLE002,PATIENT002,BATCH001,cancer_panel_v1"
        
        csv_file = tmp_path / "samples.csv"
        csv_file.write_text(csv_content)
        
        # Parse the samplesheet
        samples = parse_samplesheet(csv_file)
        
        # Verify results
        assert len(samples) == 2
        assert samples[0]["sample_id"] == "SAMPLE001"
        assert samples[0]["patient_id"] == "PATIENT001"
        assert samples[0]["batch_id"] == "BATCH001"
        assert samples[0]["panel"] == "cancer_panel_v1"
    
    def test_missing_required_header(self, tmp_path):
        """Test parsing a samplesheet with missing required headers."""
        # Create a sample CSV file with missing patient_id
        csv_content = "sample_id,batch_id,panel\n" \
                      "SAMPLE001,BATCH001,cancer_panel_v1"
        
        csv_file = tmp_path / "samples.csv"
        csv_file.write_text(csv_content)
        
        # Parsing should raise ValueError
        with pytest.raises(ValueError) as exc_info:
            parse_samplesheet(csv_file)
        
        assert "missing required headers" in str(exc_info.value)
        assert "patient_id" in str(exc_info.value)
    
    def test_file_not_found(self):
        """Test parsing a non-existent samplesheet."""
        with pytest.raises(FileNotFoundError):
            parse_samplesheet("/path/to/nonexistent/file.csv")


if __name__ == "__main__":
    unittest.main() 