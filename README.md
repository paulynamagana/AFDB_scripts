# AlphaFold Analysis Toolkit

A collection of Python scripts for retrieving, processing, and visualizing data from the AlphaFold Database, currently focused on AlphaMissense annotations with plans for expansion.

## Overview
This toolkit is being developed to enable researchers to:

Fetch protein structure data from the AlphaFold Database (AFDB) via their API
Process various protein annotations, currently focusing on AlphaMissense pathogenicity scores
Integrate annotation data with structural information from PDB files
Generate publication-quality visualizations of protein features
Analyze multiple proteins in batch from a simple input file

**Note: This project is actively being expanded to include additional data types and analysis tools beyond AlphaMissense.**

## Current Features

API Integration: Automated retrieval of data from the AlphaFold Database
AlphaMissense Analysis:

Download and process AlphaMissense pathogenicity scores
Calculate average pathogenicity scores for each residue position
Modify PDB files to incorporate AlphaMissense scores in place of B-factors


Data Visualization:

Generate protein-specific heatmaps of variant pathogenicity scores
Create plots comparing AlphaMissense scores with AlphaFold pLDDT confidence values


Batch Processing: Process multiple UniProt IDs from a single input file
Structured Output: Organize all results in a clean directory structure

## Getting Started

### Installation
Prerequisites

- Python 3.6 or higher
- pip package manager

Setup

1. Clone the repository:

```bash
git clone https://github.com/yourusername/alphafold-analysis-toolkit.git
cd alphafold-analysis-toolkit
```

2. Create and activate a virtual environment (recommended):
```bash
 python -m venv venv

# On Linux/macOS
source venv/bin/activate

# On Windows
venv\Scripts\activate
```

3. Install required dependencies:
```bash
pip install -r requirements.txt
```


### Usage
#### Step 1: Create Input File

Create a text file containing UniProt accession IDs of proteins you wish to analyze, comma-separated:
```
P04637,Q9Y2H6,P0DTD1
```

#### Step 2: Retrieve Data
First, fetch data from the AlphaFold Database API:

```bash
python alphafold_api_downloader.py uniprot_ids.txt
```


#### Step 3: Process and Visualize Data
Process the downloaded data to generate visualizations and modified PDB files:
```bash
python AM_data_processing.py uniprot_ids.txt
```

### Output Files
The toolkit generates the following outputs in the data_output directory:

CSV files: Extracted AlphaMissense data for each protein
Modified PDB files: PDB files with B-factors replaced by average pathogenicity scores
Heatmap visualizations: Pathogenicity scores for each amino acid position
Comparison plots: AlphaMissense scores vs. pLDDT confidence values

### Example Outputs

- #### AlphaMissense Pathogenicity Heatmap

<img src="https://github.com/paulynamagana/afdb-analysis-tools/blob/main/data_output/AM_heatmap_P04637.png?raw=true" alt="heatmap_P04637" style="width:500px;align:center"/>

This heatmap visualizes variant pathogenicity scores across each position in the protein sequence. The x-axis shows residue positions, with reference amino acids indicated on the top axis. The y-axis displays possible amino acid substitutions, while the color intensity represents the pathogenicity score (0-1).

- #### AlphaMissense vs. pLDDT Comparison Plot

<img src="https://github.com/paulynamagana/afdb-analysis-tools/blob/main/data_output/graph_plDDT-AM-score_P04637.png?raw=true" alt="plddt_AM_score_P04637" style="width:500px;align:center"/>

This plot compares average AlphaMissense pathogenicity scores (blue) with AlphaFold pLDDT confidence scores (orange) across the protein sequence. This comparison helps identify regions where structural uncertainty may correlate with pathogenicity potential.

- #### Google Colab Integration

For those that prefer using Google Colab:

 <a href="https://colab.research.google.com/github/paulynamagana/afdb-analysis-tools/blob/main/notebooks/alphamissense_notebook.ipynb"><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/d/d0/Google_Colaboratory_SVG_Logo.svg/1280px-Google_Colaboratory_SVG_Logo.svg.png" alt="google-colab" style="width:200px;align:center"/></a>

### Technical Implementation
#### Core Components

1. API integration layer (```alphafold_api_downloader.py```):
- Implements REST API client for AlphaFold Database
- Features robust error handling with exponential backoff retry logic
- Modular design for future extension to additional data sources

<br>

2. Data processing engine (```AM_data_processing.py```):
- Parses and transforms raw annotation data
- Implements data aggregation techniques for residue-level metrics
- Provides tools for structural file manipulation and enhancement

<br>

3. Data visualization framework (```plots.py```):
- Custom color schemes designed for optimal biological data visualization
- Dynamic scaling for different protein sequence lengths
- Dual-axis plotting for amino acid and position information

### Prerequisites

* Python 3.6 or higher

