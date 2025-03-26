#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
EHR Integration Client for Clinical Pipeline

This module provides a client for integrating genomic reports with
Electronic Health Record (EHR) systems via HL7/FHIR interfaces.
"""

import os
import sys
import json
import logging
import argparse
import datetime
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Union

import requests
from requests.auth import HTTPBasicAuth

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("EHR_CLIENT")

class EHRClient:
    """Client for integrating with Electronic Health Record systems."""
    
    def __init__(
        self, 
        api_url: str, 
        api_key: Optional[str] = None, 
        username: Optional[str] = None,
        password: Optional[str] = None,
        verify_ssl: bool = True,
        fhir_version: str = 'R4'
    ):
        """
        Initialize the EHR client.
        
        Args:
            api_url: Base URL for the EHR API
            api_key: API key for authentication (if applicable)
            username: Username for basic authentication (if applicable)
            password: Password for basic authentication (if applicable)
            verify_ssl: Whether to verify SSL certificates
            fhir_version: FHIR version to use (default: R4)
        """
        self.api_url = api_url.rstrip('/')
        self.api_key = api_key
        self.username = username
        self.password = password
        self.verify_ssl = verify_ssl
        self.fhir_version = fhir_version
        
        # Set up session
        self.session = requests.Session()
        self.session.verify = verify_ssl
        
        # Configure authentication
        if api_key:
            self.session.headers.update({'Authorization': f'Bearer {api_key}'})
        elif username and password:
            self.session.auth = HTTPBasicAuth(username, password)
        
        # Set content type
        self.session.headers.update({
            'Content-Type': 'application/fhir+json',
            'Accept': 'application/fhir+json'
        })
    
    def send_diagnostic_report(
        self, 
        report_data: Dict, 
        patient_id: str, 
        report_pdf: Optional[Union[str, Path]] = None
    ) -> Dict:
        """
        Send a diagnostic report to the EHR system.
        
        Args:
            report_data: Dictionary containing report data
            patient_id: Patient identifier in the EHR system
            report_pdf: Path to PDF report file (optional)
            
        Returns:
            Dictionary containing the response from the EHR system
        """
        # Create FHIR DiagnosticReport resource
        fhir_report = self._create_fhir_diagnostic_report(report_data, patient_id)
        
        # Send report to FHIR endpoint
        response = self.session.post(
            f"{self.api_url}/DiagnosticReport",
            json=fhir_report
        )
        
        # Check response
        response.raise_for_status()
        result = response.json()
        
        # If report PDF is provided, attach it as a DocumentReference
        if report_pdf:
            doc_ref = self._attach_pdf_document(result['id'], patient_id, report_pdf)
            result['documentReference'] = doc_ref
        
        return result
    
    def _create_fhir_diagnostic_report(
        self, 
        report_data: Dict, 
        patient_id: str
    ) -> Dict:
        """
        Create a FHIR DiagnosticReport resource from report data.
        
        Args:
            report_data: Dictionary containing report data
            patient_id: Patient identifier in the EHR system
            
        Returns:
            Dictionary containing the FHIR DiagnosticReport resource
        """
        # Basic report structure
        report = {
            "resourceType": "DiagnosticReport",
            "id": f"genomic-report-{uuid.uuid4()}",
            "status": "final",
            "category": [
                {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                            "code": "GE",
                            "display": "Genetics"
                        }
                    ]
                }
            ],
            "code": {
                "coding": [
                    {
                        "system": "http://loinc.org",
                        "code": "81247-9",
                        "display": "Master HL7 genetic variant reporting panel"
                    }
                ],
                "text": report_data.get('report_title', 'Genomic Variant Report')
            },
            "subject": {
                "reference": f"Patient/{patient_id}"
            },
            "effectiveDateTime": report_data.get('report_date', datetime.datetime.now().isoformat()),
            "issued": datetime.datetime.now().isoformat(),
            "performer": [
                {
                    "display": f"{report_data.get('reviewer_name', 'Unknown')}, {report_data.get('reviewer_credentials', '')}"
                }
            ],
            "conclusion": report_data.get('clinical_interpretation', 'No interpretation provided.')
        }
        
        # Add observations for each variant
        if 'variants' in report_data:
            report["result"] = []
            
            for variant in report_data['variants']:
                observation_id = f"variant-{uuid.uuid4()}"
                
                # Create Observation resource for variant
                observation = {
                    "resourceType": "Observation",
                    "id": observation_id,
                    "status": "final",
                    "category": [
                        {
                            "coding": [
                                {
                                    "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                    "code": "laboratory",
                                    "display": "Laboratory"
                                }
                            ]
                        }
                    ],
                    "code": {
                        "coding": [
                            {
                                "system": "http://loinc.org",
                                "code": "69548-6",
                                "display": "Genetic variant assessment"
                            }
                        ]
                    },
                    "subject": {
                        "reference": f"Patient/{patient_id}"
                    },
                    "effectiveDateTime": report_data.get('report_date', datetime.datetime.now().isoformat()),
                    "valueCodeableConcept": {
                        "coding": [
                            {
                                "system": "http://loinc.org",
                                "code": self._map_classification_to_loinc(variant.get('classification', 'Uncertain Significance')),
                                "display": variant.get('classification', 'Uncertain Significance')
                            }
                        ]
                    },
                    "component": [
                        {
                            "code": {
                                "text": "Gene studied"
                            },
                            "valueString": variant.get('gene', '')
                        },
                        {
                            "code": {
                                "text": "Variant HGVS name"
                            },
                            "valueString": variant.get('hgvs', '')
                        },
                        {
                            "code": {
                                "text": "Zygosity"
                            },
                            "valueString": variant.get('zygosity', '')
                        },
                        {
                            "code": {
                                "text": "Clinical significance"
                            },
                            "valueString": variant.get('significance', '')
                        }
                    ]
                }
                
                # Add reference to the observation
                report["result"].append({
                    "reference": f"#{observation_id}"
                })
                
                # Include the observation as a contained resource
                if not "contained" in report:
                    report["contained"] = []
                report["contained"].append(observation)
        
        return report
    
    def _attach_pdf_document(
        self, 
        report_id: str, 
        patient_id: str, 
        pdf_path: Union[str, Path]
    ) -> Dict:
        """
        Attach a PDF document to a report in the EHR system.
        
        Args:
            report_id: Report identifier in the EHR system
            patient_id: Patient identifier in the EHR system
            pdf_path: Path to PDF file
            
        Returns:
            Dictionary containing the response from the EHR system
        """
        # Encode PDF as base64
        with open(pdf_path, 'rb') as pdf_file:
            import base64
            pdf_data = base64.b64encode(pdf_file.read()).decode('utf-8')
        
        # Create DocumentReference resource
        doc_ref = {
            "resourceType": "DocumentReference",
            "status": "current",
            "docStatus": "final",
            "type": {
                "coding": [
                    {
                        "system": "http://loinc.org",
                        "code": "81247-9",
                        "display": "Genetic testing report"
                    }
                ]
            },
            "subject": {
                "reference": f"Patient/{patient_id}"
            },
            "date": datetime.datetime.now().isoformat(),
            "description": "Genomic Variant Report (PDF)",
            "content": [
                {
                    "attachment": {
                        "contentType": "application/pdf",
                        "data": pdf_data,
                        "title": f"Genomic_Report_{report_id}.pdf"
                    }
                }
            ],
            "context": {
                "related": [
                    {
                        "reference": f"DiagnosticReport/{report_id}"
                    }
                ]
            }
        }
        
        # Send document to FHIR endpoint
        response = self.session.post(
            f"{self.api_url}/DocumentReference",
            json=doc_ref
        )
        
        # Check response
        response.raise_for_status()
        return response.json()
    
    def _map_classification_to_loinc(self, classification: str) -> str:
        """Map variant classification to LOINC code."""
        mapping = {
            'Pathogenic': 'LA6668-3',
            'Likely Pathogenic': 'LA26332-9',
            'Uncertain Significance': 'LA26333-7',
            'Likely Benign': 'LA26334-5',
            'Benign': 'LA6675-8'
        }
        
        return mapping.get(classification, 'LA26333-7')  # Default to VUS
    
    def get_patient_demographics(self, patient_id: str) -> Dict:
        """
        Get patient demographics from the EHR system.
        
        Args:
            patient_id: Patient identifier in the EHR system
            
        Returns:
            Dictionary containing patient demographics
        """
        response = self.session.get(f"{self.api_url}/Patient/{patient_id}")
        response.raise_for_status()
        return response.json()
    
    def query_patient_by_mrn(self, medical_record_number: str) -> Dict:
        """
        Query for a patient by medical record number.
        
        Args:
            medical_record_number: Patient's medical record number
            
        Returns:
            Dictionary containing patient information
        """
        params = {
            'identifier': f"MRN|{medical_record_number}"
        }
        
        response = self.session.get(f"{self.api_url}/Patient", params=params)
        response.raise_for_status()
        
        result = response.json()
        
        # Check if patient was found
        if result.get('total', 0) == 0:
            raise ValueError(f"No patient found with MRN {medical_record_number}")
        
        return result['entry'][0]['resource']
    
    def verify_patient_exists(self, patient_id: str) -> bool:
        """
        Verify if a patient exists in the EHR system.
        
        Args:
            patient_id: Patient identifier in the EHR system
            
        Returns:
            True if patient exists, False otherwise
        """
        try:
            response = self.session.get(f"{self.api_url}/Patient/{patient_id}")
            response.raise_for_status()
            return True
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                return False
            raise


def send_report_to_ehr(
    report_json_path: Union[str, Path],
    patient_id: str,
    report_pdf_path: Optional[Union[str, Path]] = None,
    ehr_api_url: str = 'https://fhir-api.hospital.org/fhir',
    api_key: Optional[str] = None,
    username: Optional[str] = None,
    password: Optional[str] = None
) -> None:
    """
    Send a clinical report to an EHR system.
    
    Args:
        report_json_path: Path to report JSON file
        patient_id: Patient identifier in the EHR system
        report_pdf_path: Path to report PDF file (optional)
        ehr_api_url: URL of the EHR FHIR API
        api_key: API key for authentication (if applicable)
        username: Username for basic authentication (if applicable)
        password: Password for basic authentication (if applicable)
    """
    logger.info(f"Sending report to EHR for patient {patient_id}")
    
    # Load report data
    with open(report_json_path, 'r') as f:
        report_data = json.load(f)
    
    # Create EHR client
    client = EHRClient(
        api_url=ehr_api_url,
        api_key=api_key,
        username=username,
        password=password
    )
    
    # Verify patient exists
    if not client.verify_patient_exists(patient_id):
        raise ValueError(f"Patient {patient_id} not found in EHR system")
    
    # Send report
    result = client.send_diagnostic_report(
        report_data=report_data,
        patient_id=patient_id,
        report_pdf=report_pdf_path
    )
    
    logger.info(f"Report successfully sent to EHR. Report ID: {result.get('id')}")
    
    # Return the result
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Send genomic report to EHR")
    parser.add_argument('--report', required=True, help='Path to report JSON file')
    parser.add_argument('--patient_id', required=True, help='Patient ID in the EHR system')
    parser.add_argument('--pdf', help='Path to report PDF file (optional)')
    parser.add_argument('--api_url', default=os.environ.get('EHR_API_URL'), help='EHR FHIR API URL')
    parser.add_argument('--api_key', default=os.environ.get('EHR_API_KEY'), help='API key for authentication')
    parser.add_argument('--username', default=os.environ.get('EHR_USERNAME'), help='Username for basic authentication')
    parser.add_argument('--password', default=os.environ.get('EHR_PASSWORD'), help='Password for basic authentication')
    
    args = parser.parse_args()
    
    # Check API URL
    if not args.api_url:
        parser.error('EHR API URL is required (provide via --api_url or EHR_API_URL environment variable)')
    
    # Check authentication
    if not args.api_key and not (args.username and args.password):
        parser.error('Authentication required (provide API key or username/password)')
    
    try:
        send_report_to_ehr(
            report_json_path=args.report,
            patient_id=args.patient_id,
            report_pdf_path=args.pdf,
            ehr_api_url=args.api_url,
            api_key=args.api_key,
            username=args.username,
            password=args.password
        )
        
    except Exception as e:
        logger.error(f"Error sending report to EHR: {str(e)}")
        sys.exit(1) 