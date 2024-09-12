# data_processing.py
import pandas as pd
import numpy as np
import requests
import os
import logging
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
### Change the logging to DEBUG when needed to turn the logs back on

def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        logging.info(f"Created directory: {directory}")


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

    
def plot_am_heatmap(am_data):
    """
    Plots a heatmap for the AlphaMissense data with reference residues on x and alternative_aa on y axis
    Colour coded following the original colours from the AlphaFold Database, see: https://alphafold.ebi.ac.uk/entry/Q5VSL9
    """

    #custom palette
    def create_custom_colormap():
        cdict = {
            'red': [
                (0.0, 56/255, 56/255),
                (0.34, 204/255, 204/255),
                (0.464, 204/255, 204/255),
                (1.0, 165/255, 165/255)
            ],
            'green': [
                (0.0, 83/255, 83/255),
                (0.34, 204/255, 204/255),
                (0.464, 204/255, 204/255),
                (1.0, 13/255, 13/255)
            ],
            'blue': [
                (0.0, 163/255, 163/255),
                (0.34, 204/255, 204/255),
                (0.464, 204/255, 204/255),
                (1.0, 18/255, 18/255)
            ]
        }
        return LinearSegmentedColormap('CustomMap', segmentdata=cdict)

    # custom colormap
    custom_cmap = create_custom_colormap()

    # pivot table
    pivot_table = am_data.pivot_table(values='pathogenicity_score', index='alternative_aa', columns='reference_aa', aggfunc='mean')
    pivot_table = pd.pivot_table(am_data, values='pathogenicity_score',
    index='alternative_aa', columns='residue_number')

    #custom_order = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "P", "A", "V", "I", "L", "M", "G", "F","Y","W"]

    # Reindex the pivot table
    #pivot_table = pivot_table.reindex(custom_order)

    plt.figure(figsize=(20, 6))

    ax = sns.heatmap(pivot_table, cmap=custom_cmap, vmin=0, vmax=1,
    cbar_kws={'label': 'AlphaMissense score'}) # Limits for the color scale

    ax.set_xlabel('Residue Number')
    ax.set_ylabel('Alternative Amino Acid')
    plt.title('AlphaMissense Pathogenicity Heatmap ')

    xticks = range(0, pivot_table.shape[1], 50)
    ax.set_xticks(xticks)
    ax.set_xticklabels(pivot_table.columns[xticks])
    ax.set_facecolor('black') #Set background black for matching AA
    plt.yticks(rotation = 0)

    cbar = ax.collections[0].colorbar
    cbar.set_ticks([i / 10.0 for i in range(11)])
    cbar.set_ticklabels([f'{i / 10.0:.1f}' for i in range(11)])
    plt.tight_layout()


    # Save the plot to a file
    output_directory = "data_output"
    ensure_directory_exists(output_directory)
    
    output_file = os.path.join(output_directory, f"AM_heatmap.png")
    plt.savefig(output_file) #save file 
    plt.show() # Show plot

    plt.close() # Close the figure after saving to free up resources



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
