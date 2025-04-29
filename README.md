# AlphaFold Database Scripts

This repository contains two Python scripts for working with the AlphaFold Database API:

* `alphafold_api_downloader.py`
* `AM_data_processing.py`

These scripts can be used to download AlphaMissense (AM) data and PDB URLs for a list of UniProt IDs, calculate average pathogenicity scores for each residue, modify PDB files with these scores, and visualize the data using heatmaps.

## Getting Started

### Prerequisites

* Python 3.6 or higher
* A `requirements.txt` file is provided to easily create a dedicated environment with the necessary Python libraries.

### Installation

1.  Clone this repository (or download the script file and the `requirements.txt` file).
2.  Navigate to the repository directory in your terminal or command prompt.
3.  Create a virtual environment (recommended) using your preferred method (e.g., `venv`, `conda`). For example, using `venv`:
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Linux/macOS
    venv\Scripts\activate  # On Windows
    ```
4.  Install the required Python libraries from the `requirements.txt` file:
    ```bash
    pip install -r requirements.txt
    ```


## `alphafold_api_downloader.py`

### Usage
* Prepare your input file:

Create a text file (e.g., uniprot_ids.txt) containing the UniProt IDs you want to process, one per line.
You can separate the IDs with commas if you prefer.

* Run the script:
`python alphafold_api_downloader.py <filename> [options]`

Replace <filename> with the name of your input file.

* Options:

--extract_am_files: Extract files from amAnnotationsHg19Url and amAnnotationsHg38Url (default: False).

Example:

`python alphafold_api_downloader.py uniprot_ids.txt --extract_am_files`

### Output

The script will create a directory named downloaded_files in the same directory as the script.
The downloaded files will be saved in this directory.
By default, the script will only download the AlphaMissense URL for each UniProt ID.
If you use the `--extract_am_files` option, it will download all files.

## `AM_data_processing.py`

The script performs the following steps:

1.  **Fetches Data:** Retrieves AlphaMissense annotation URLs and PDB URLs from the AlphaFold Database API using a provided UniProt ID.
2.  **Extracts AlphaMissense Data:** Downloads and parses the AlphaMissense data, extracting pathogenicity scores and corresponding residue information.
3.  **Calculates Average Pathogenicity:** Computes the average pathogenicity score for each residue across all observed variants.
4.  **Modifies PDB File:** Downloads the PDB file and writes a new PDB file where the B-factor column (columns 61-66) for each ATOM and HETATM record is replaced with the calculated average AlphaMissense pathogenicity score for that residue.
5.  **Generates Plots:**
    * Creates a heatmap visualizing the pathogenicity scores for each variant across the protein sequence.
    * Generates a scatter plot comparing the average pathogenicity scores (stored in the modified PDB) against the pLDDT values extracted from the original PDB file.

### Usage
### Running the Script

1.  Create a text file (e.g., `uniprot_ids.txt`) containing one or more UniProt IDs, separated by commas (no spaces). For example:
    ```
    P0DTD1,Q9Y2H6
    ```
2.  Ensure your virtual environment is activated.
3.  Run the script, providing the text file as an argument:
    ```bash
    python your_script_name.py uniprot_ids.txt
    ```
    Replace `your_script_name.py` with the actual name of the Python script.


### Output

Upon successful execution, the script will create a directory named `data_output` containing:

* CSV files of the downloaded AlphaMissense data.
* Modified PDB files named `AM_scores_AF-{UniProtID}-F1-model_v4.pdb`, where the B-factor column contains the average AlphaMissense pathogenicity scores.
* A heatmap plot (PNG file) visualizing the pathogenicity scores for each variant.
* A scatter plot (PNG file) comparing the average pathogenicity scores and pLDDT values.

![AM heatmap](https://github.com/paulynamagana/AFDB_scripts/blob/main/data_output/AM_heatmap_P04637.png)

![pLDDT and AM scores](https://github.com/paulynamagana/AFDB_scripts/blob/main/data_output/graph_plDDT-AM-score_P04637.png)


## Accessibility via Google Colab

You can also run this script easily on Google Colab:

[![Open In Colab](https://colab.research.google.com/github/paulynamagana/afdb-analysis-tools/blob/main/notebooks/alphamissense_notebook.ipynb)
