# BLAST Protein Search & 3D Visualization

## Overview

This is a Streamlit-based bioinformatics application that enables researchers to search for similar protein structures using BLAST (Basic Local Alignment Search Tool) and visualize the results in 3D. The application integrates with the NCBI BLAST service to search against the Protein Data Bank (PDB) database, downloads matching protein structures, and provides interactive 3D visualization capabilities. The tool is designed for protein analysis workflows, allowing users to input protein sequences, perform similarity searches, and explore structural homologs through an intuitive web interface.

## Recent Changes

### 2025-09-26: Project Import and Setup
- Successfully imported GitHub project into Replit environment
- Configured Streamlit server settings for Replit (0.0.0.0:5000 with headless mode)
- Set up workflow for continuous development server
- Configured deployment settings for autoscale production deployment
- All dependencies properly installed and working
- Application running successfully on port 5000

## User Preferences

Preferred communication style: Simple, everyday language.

## System Architecture

### Frontend Architecture
The application uses Streamlit as the primary web framework, providing a Python-based approach to building interactive web applications. The interface is organized into multiple tabs for different workflow stages: sequence input, BLAST results, 3D structures, and analysis/export. The sidebar contains configurable search parameters including E-value thresholds and hit limits. Session state management is implemented to maintain data persistence across user interactions and tab navigation.

### Backend Architecture
The backend follows a modular utility-based architecture with three main components:
- **BLAST utilities** (`blast_utils.py`): Handles protein sequence searches against NCBI databases using BioPython's BLAST interface
- **PDB utilities** (`pdb_utils.py`): Manages downloading and processing of 3D protein structures from the RCSB PDB database
- **Visualization utilities** (`visualization.py`): Creates interactive 3D molecular visualizations using py3Dmol and Plotly

The application uses BioPython for biological data processing and sequence handling, providing standardized methods for working with protein sequences and structure files.

### Data Processing Pipeline
The workflow follows a sequential pipeline: sequence input → BLAST search → PDB ID extraction → structure download → 3D visualization. Results are cached in session state to prevent redundant API calls. The application supports both manual sequence input and file uploads for batch processing.

### Visualization Strategy
3D molecular visualization is implemented using py3Dmol for interactive structure viewing with cartoon and surface representations. Plotly is used for additional data visualization and alignment plots. The visualization components are embedded directly in the Streamlit interface for seamless user experience.

## External Dependencies

### Bioinformatics Services
- **NCBI BLAST API**: Remote protein sequence similarity searching via BioPython's NCBIWWW interface
- **RCSB Protein Data Bank**: Source for downloading 3D protein structure files in PDB format

### Python Libraries
- **Streamlit**: Web application framework and user interface
- **BioPython**: Biological sequence analysis, BLAST integration, and PDB file parsing
- **py3Dmol**: 3D molecular structure visualization in web browsers
- **Plotly**: Interactive plotting and data visualization
- **Pandas/NumPy**: Data manipulation and numerical computing

### File System
- Local storage for downloaded PDB structure files in `pdb_files/` directory
- Temporary file handling for user uploads and processed sequences

### API Integrations
The application makes HTTP requests to external bioinformatics databases and services, with error handling for network timeouts and service availability. Multiple PDB download sources are configured as fallbacks to ensure structure retrieval reliability.