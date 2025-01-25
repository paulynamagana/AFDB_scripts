
import pandas as pd
import numpy as np
import requests
import os
import logging
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from alphafold_api_downloader import fetch_AFDB_data
import plots
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
### Change the logging to DEBUG when needed to turn the logs back on

def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        logging.info(f"Created directory: {directory}")


def extract_alpha_missense_url(data):
    """Extracts the AlphaMissense URL from the API data.

    Args:
        data (dict): The JSON data from the API.

    Returns:
        str: The AlphaMissense URL, or an error message.
    """
    if data:
        logging.info(f"Retrieving AM data")
        return data[0].get('amAnnotationsUrl', f"No AlphaMissense data for this protein.")
    else:
        return "Error: No data provided."


def extract_pdb_url(data):
    """Extracts the PDB URL from the API data.

    Args:
        data (dict): The JSON data from the API.

    Returns:
        str: The PDB url, or an error message.
    """
    if data:
        logging.info(f"Retrieving url for PDB file")
        return data[0].get('pdbUrl', "Failed to retrieve PDB URL.")
    else:
        return "Error: No data provided."
    


def extract_am_data(am_url):
    """
    Extract AlphaMissense data from the url and saves it as a CSV file
    """
    try:
        am_file = pd.read_csv(am_url)
    except Exception as e:
        logging.error(f"Error reading AM file: {e}")
        return None
        
    reference_aa = am_file["protein_variant"].str.extract(r'^([A-Z])')[0]
    alternative_aa = am_file["protein_variant"].str.extract('([A-Z])$')[0]
    residue_number = pd.to_numeric(am_file["protein_variant"].str.extract(r'([0-9]+)')[0])
    pathogenicity_score = pd.to_numeric(am_file['am_pathogenicity'])
        
    am_data = pd.DataFrame({
        'reference_aa': reference_aa,
        'residue_number': residue_number,
        "alternative_aa":alternative_aa,
        "pathogenicity_score": pathogenicity_score})
        
    output_directory = "data_output"
    ensure_directory_exists(output_directory)
    output_filename = os.path.basename(am_url)
    output_file = os.path.join(output_directory, output_filename)

    
    try:
        am_data.to_csv(output_file, index=False)
        logging.info(f"AM data saved to: {output_filename}")
    except Exception as e:
        logging.error(f"Error saving AM data: {e}")
        
    return am_data
        
def calculate_average_pathogenicity(am_data):
    """Calculates average pathogenicity scores per residue from AM data
    """
    if am_data is None:
        logging.error("No AlphaMissense data available")
        return None
    
    grouped = am_data.groupby(['residue_number'])['pathogenicity_score'].mean().reset_index()
    max_residue_number = grouped['residue_number'].max()
    average_scores = np.full(max_residue_number + 1, np.nan)
    
    for _, row in grouped.iterrows():
        residue_number = int(row['residue_number'])
        average_scores[residue_number] = round(row['pathogenicity_score'], 4)
    return average_scores


def modify_pdb_with_am_data(pdb_url, average_scores_file):
    
    if pdb_url:
        try:
            logging.info("Retrieving PDB file")
            pdb_content = requests.get(pdb_url).text.splitlines()
        except requests.exceptions.RequestException as e:
            logging.error(f"Failed to retrieve the PDB file: {e}")  

    output_directory = "data_output"
    ensure_directory_exists(output_directory)
    
    base_filename = os.path.basename(pdb_url)
    output_filename = f"AM_scores_{base_filename}"
    output_file = os.path.join(output_directory, output_filename)
    
    try:
        logging.info(f"Writing modified PDB data to: {output_file}")
        with open(output_file, "w", encoding= "utf-8") as out_file:
            for line in pdb_content:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    residue_number = int(line[22:26].strip())
                    if residue_number < len(average_scores_file) and not np.isnan(average_scores_file[residue_number]):
                        value = average_scores_file[residue_number]
                        value_str = f"{value:.2f}"
                        while len(value_str) < 6:
                            value_str = " " + value_str
                        edit_line = line[:60] + value_str + line[66:]
                        out_file.write(edit_line + '\n')
                    else:
                        out_file.write(line + '\n')
                else:
                    out_file.write(line + '\n')
    except IOError as e:
        logging.error(f"Error writing to file: {e}")


def extract_pathogenicity_and_plddt(am_file_path, pdb_file_url):
    """
    Extracts pathogenicity scores from an AM file and pLDDT values from a PDB file URL.

    Args:
    am_file_path (str): Path to the AM file containing pathogenicity scores.
    pdb_file_url (str): URL of the PDB file containing pLDDT values.

    Returns:
    tuple: Two lists:
    - pathogenicity_scores (list): List of pathogenicity scores.
    - plddt_scores (list): List of pLDDT scores.
    """
    pathogenicity_scores = []
    plddt_scores = []

    # Extract Pathogenicity Scores from AM File
    try:
        with open(am_file_path, "r") as f:
            am_pdb = f.read()

        for line in am_pdb.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                pathogenicity_scores.append(float(line[60:66].strip()))

    except FileNotFoundError:
        logging.error(f"AM file not found: {am_file_path}")
    except ValueError:
        logging.error("Error parsing pathogenicity scores in the AM file.")

    # Extract pLDDT Scores from PDB URL
    if pdb_file_url:
        logging.info(f"Fetching PDB data from {pdb_file_url}")

        try:
            response = requests.get(pdb_file_url)
            response.raise_for_status()
            pdb_content = response.text.splitlines()

            for line in pdb_content:
                if line.startswith(("ATOM", "HETATM")):
                    plddt_scores.append(float(line[60:66].strip()))
        
        except requests.exceptions.RequestException as e:
            logging.error(f"Error fetching PDB data: {e}")
        except ValueError:
            logging.error("Error parsing pLDDT scores in the PDB file.")

    return pathogenicity_scores, plddt_scores


def process_am_from_file(filename):
    """
    Reads UniProt IDs from a text file, fetches data and extracts URLs.

    Args:
        filename (str): The name of the file containing UniProt IDs.
        extract_am_files (bool): Whether to extract files from amAnnotationsHg19Url and amAnnotationsHg38Url.
    """
    
    with open(filename, 'r') as f:
        uniprot_ids = f.read().replace(' ', '').split(',')  # Remove whitespace and split by comma
        
    output_dir = "downloaded_files"
    os.makedirs(output_dir, exist_ok=True)

    
    for uniprot_id in uniprot_ids:
        try:
            data = fetch_AFDB_data(uniprot_id)
            if data:
                am_data_url = extract_alpha_missense_url(data)
                
                if not am_data_url or not am_data_url.startswith("http"):  # Check for both None and invalid URLs
                    continue

                am_data = extract_am_data(am_data_url)
                pdb_data_url = extract_pdb_url(data)
                
                if am_data is not None:
                    average_scores = calculate_average_pathogenicity(am_data)
                    modify_pdb_with_am_data(pdb_data_url, average_scores)
                    print(uniprot_id)
                    plots.plot_am_heatmap(am_data, uniprot_id)
                

                file_path= f"data_output/AM_scores_AF-{uniprot_id}-F1-model_v4.pdb"
                pathogenicity_scores, plddt_scores = extract_pathogenicity_and_plddt(file_path, pdb_data_url)

                print(uniprot_id)
                plots.plot_scores(pathogenicity_scores, plddt_scores, uniprot_id)
                            
        except requests.exceptions.RequestException as e:
            print(f"Error processing {uniprot_id}: {e}")



# Define the main function
def main():
    parser = argparse.ArgumentParser(description='Fetch and extract URLs from the AlphaFold Database API.')
    parser.add_argument('filename', help='Text file containing UniProt IDs, one per line')
    args = parser.parse_args()

    # Call the process_uniprot_ids_from_file function with the provided filename and output option
    process_am_from_file(args.filename)

# Call the main function if this script is run directly
if __name__ == '__main__':
    main()