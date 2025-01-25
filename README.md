AlphaFold Database Scripts
This repository contains two Python scripts for working with the AlphaFold Database API:

alphafold_api_downloader.py
AM_data_processing.py
These scripts can be used to download AlphaMissense (AM) data and PDB URLs for a list of UniProt IDs, calculate average pathogenicity scores for each residue, modify PDB files with these scores, and visualize the data using heatmaps.

Prerequisites
Python 3.6 or higher
requests library
Installation
Clone the repository: git clone [invalid URL removed]
Install the required library: pip install requests
Usage
alphafold_api_downloader.py
Prepare your input file:

Create a text file (e.g., uniprot_ids.txt) containing the UniProt IDs you want to process, one per line.
You can separate the IDs with commas if you prefer.

Run the script:
python alphafold_api_downloader.py <filename> [options]

Replace <filename> with the name of your input file.

Options:

--extract_am_files: Extract files from amAnnotationsHg19Url and amAnnotationsHg38Url (default: False).

Example:

python alphafold_api_downloader.py uniprot_ids.txt --extract_am_files

Output

The script will create a directory named downloaded_files in the same directory as the script.
The downloaded files will be saved in this directory.
By default, the script will only download the AlphaMissense URL for each UniProt ID.
If you use the --extract_am_files option, it will download all files.

AM_data_processing.py
This script processes UniProt IDs to extract AlphaMissense (AM) data and PDB URLs from the AlphaFold Database API. It then calculates average pathogenicity scores for each residue, modifies PDB files with these scores, and plots heatmaps to visualize the data.

Usage

To use this script, you need to provide a text file containing a list of UniProt IDs. You can then run the script using the following command:

python AM_data_processing.py <filename>

where <filename> is the name of the text file containing the UniProt IDs.

Output

The script will save the AM data and modified PDB files to the data_output directory. It will also display the UniProt ID and plot the AM heatmap and score plots for each UniProt ID processed.

[![DOI](https://zenodo.org/badge/834257572.svg)](https://doi.org/10.5281/zenodo.13942608)
