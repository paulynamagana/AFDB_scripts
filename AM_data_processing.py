# data_processing.py

import pandas as pd
import numpy as np
import requests
import os
import logging
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

from api_handler import retrieve_AFDB

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
#Change the logging to DEBUG when needed to turn the logs back on

def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        logging.info(f"Created directory: {directory}")


class extract_AM_data(retrieve_AFDB):
    """
    Extract AlphaMissense data from the previous call and store it as a df with 3 columns:
    reference_aa, residue_number and pathogenicity_score
    """
    def __init__(self, AM_url, uniprot_accession):
        self.AM_url = AM_url
        retrieve_AFDB.__init__(self, uniprot_accession)
        self.uniprot_accession = uniprot_accession
        
    def extract_AM(self):
        try:
            AM_file = pd.read_csv(self.AM_url)
        except Exception as e:
            logging.error(f"Error reading AM file: {e}")
            return None
        
        reference_aa = AM_file["protein_variant"].str.extract(r'^([A-Z])')[0]
        alternative_aa = AM_file["protein_variant"].str.extract('([A-Z])$')[0]
        residue_number = pd.to_numeric(AM_file["protein_variant"].str.extract(r'([0-9]+)')[0])
        pathogenicity_score = pd.to_numeric(AM_file['am_pathogenicity'])
        
        AM_data = pd.DataFrame({
                'reference_aa': reference_aa,
                'residue_number': residue_number,
                "alternative_aa":alternative_aa,
                "pathogenicity_score": pathogenicity_score})
        
        output_directory = "data_output"
        ensure_directory_exists(output_directory)
        output_file = os.path.join(output_directory, f"{self.uniprot_accession}_AM_scores.csv")
      
        try:
            AM_data.to_csv(output_file, index=False)
            logging.info(f"AM data saved to: {output_file}")
        except Exception as e:
            logging.error(f"Error saving AM data: {e}")
        
        return AM_data
        


class process_AM_data(extract_AM_data, retrieve_AFDB):
    """
    Class for processing AlphaMissense data such as extracting the data, modifying the PDB file and plotting the heatmap
    """
    def __init__(self, AM_url, uniprot_accession):
        extract_AM_data.__init__(self, AM_url, uniprot_accession)
        retrieve_AFDB.__init__(self, uniprot_accession)
        
        self.AM_data_file = self.extract_AM()
        self.uniprot_accession = uniprot_accession
        
    def calculate_average_pathogenicity(self):
        if self.AM_data_file is None:
            logging.error("No AM data available")
            return None
            
        grouped = self.AM_data_file.groupby(['residue_number'])['pathogenicity_score'].mean().reset_index()
        max_residue_number = grouped['residue_number'].max()
        average_scores = np.full(max_residue_number + 1, np.nan)
        for _, row in grouped.iterrows():
            residue_number = int(row['residue_number'])
            average_scores[residue_number] = round(row['pathogenicity_score'], 4)
        return average_scores
        
    def modify_pdb_with_pathogenicity(self, pdb_file_url):
        average_scores_file = self.calculate_average_pathogenicity()
        pdb_url = self.get_pdb_url()

        if pdb_url:
            logging.info(f"Fetching PDB data from {pdb_url}")
            
            try:
                pdb_content = requests.get(pdb_url).text.splitlines()
            except requests.exceptions.RequestException as e:
                logging.error(f"Error fetching PDB data: {e}")
                return   
            
            output_directory = "data_output"
            ensure_directory_exists(output_directory)
            
            output_file = os.path.join(output_directory,f"{self.uniprot_accession}_with_AM_scores.PDB")
            logging.info(f"Writing modified PDB data to: {output_file}")
            
            try:
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
        
        else:
            logging.error("Failed to retrieve PDB URL.")



class plotting_AM_heatmap(extract_AM_data, retrieve_AFDB):
    """
    Class for plotting a heatmap for the AlphaMissense data with reference residues on x and alternative_aa on y axis
    Colour coded following the original colours from the AlphaFold Database, see: https://test.alphafold.ebi.ac.uk/entry/Q5VSL9
    """
    def __init__(self, AM_url, uniprot_accession):
        extract_AM_data.__init__(self, AM_url, uniprot_accession)
        retrieve_AFDB.__init__(self, uniprot_accession)
        
        self.AM_data_file = self.extract_AM()
        self.uniprot_accession = uniprot_accession
    
    def plot_AM_heatmap(self):
        """This code will plot the data for the UniProt accession provided and will downaload it as a png file to the local computer
        """
        pivot_table = self.AM_data_file.pivot_table(values='pathogenicity_score', index='alternative_aa', columns='reference_aa', aggfunc='mean')

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
        pivot_table = pd.pivot_table(self.AM_data_file, values='pathogenicity_score',
                                    index='alternative_aa', columns='residue_number')

        #custom_order = ["R", "H", "K", "D", "E", "S", "T", "N", "Q",  "C", "P", "A", "V", "I", "L", "M", "G", "F","Y","W"]

        # Reindex the pivot table
        #pivot_table = pivot_table.reindex(custom_order)

        plt.figure(figsize=(20, 6))

        ax = sns.heatmap(pivot_table, cmap=custom_cmap, vmin=0, vmax=1,
                         cbar_kws={'label': 'AlphaMissense score'})  # Limits for the color scale

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
        
        output_file = os.path.join(output_directory, f"pathogenicity_heatmap_{self.uniprot_accession}_AM_heatmap.png")
        plt.savefig(output_file)
        plt.close()