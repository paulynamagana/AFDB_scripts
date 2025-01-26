import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from matplotlib.colors import LinearSegmentedColormap
import logging

def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        logging.info(f"Created directory: {directory}")


def plot_plddt_legend(dpi=200):
    """
    Creates and returns a Matplotlib figure containing a legend for pLDDT (predicted Local Distance Difference Test) scores.

    Args:
        dpi (int, optional): The resolution of the figure in dots per inch. Defaults to 100.

    Returns:
        matplotlib.figure.Figure: The figure object containing the legend.
    """
    plddt_labels = ["plDDT:", "Very low (<50)", "Low (60)", "Confident (80)", "Very high (>90)"]
    plddt_colors = ["#FFFFFF", "#FF7D45", "#FFFF00", "#65CBF3", "#0000FF"]

    fig, ax = plt.subplots(figsize=(1, 0.1), dpi=dpi)  # Create figure and axis objects

    # Create legend handles (dummy bars)
    for color in plddt_colors:
        ax.bar(0, 0, color=color)

    legend = ax.legend(
        plddt_labels,
        frameon=False,
        loc="center",
        ncol=6,
        handletextpad=1,
        columnspacing=1,
        markerscale=0.5,
    )

    ax.axis("off")  # Turn off axis labels and ticks

    return fig


def plot_scores(pathogenicity_scores, plddt_scores, uniprot_id):
    """
    Plots pathogenicity and rescaled pLDDT scores against residue number using Seaborn.

    Args:
        pathogenicity_scores (list): List of pathogenicity scores.
        plddt_scores (list): List of pLDDT scores.
    """

    # Rescale pLDDT scores from 0-100 to 0-1
    plddt_rescaled = [score / 100 for score in plddt_scores]

    # Create a DataFrame for easy plotting with Seaborn
    data = {
        'Residue number': range(1, len(pathogenicity_scores) + 1),
        'Average pathogenicity': pathogenicity_scores,
        'pLDDT (rescaled)': plddt_rescaled
    }
    df = pd.DataFrame(data)

    # Create line plots with Seaborn
    plt.figure(figsize=(12, 5))
    sns.lineplot(x='Residue number', y='Average pathogenicity', data=df, label='Average pathogenicity')
    sns.lineplot(x='Residue number', y='pLDDT (rescaled)', data=df, label='pLDDT (rescaled)')

    # Add labels, title, and legend
    plt.xlabel('Residue number')
    plt.ylabel('Score')
    plt.title(f'Average AM and pLDDT scores per position ({uniprot_id})')
    plt.legend()
    plt.grid(axis='y', linestyle='--')  # Optional grid

        # Save the plot to a file
    output_directory = "data_output"
    output_file = os.path.join(output_directory, f"graph_plDDT-AM-score_{uniprot_id}.png")
    plt.savefig(output_file) #save file 
    
    # Show the plot
    #plt.show()
    plt.close()


def plot_am_heatmap(am_data, uniprot_id):
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
    plt.title(f'AlphaMissense Pathogenicity Heatmap ({uniprot_id})')

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
    
    output_file = os.path.join(output_directory, f"AM_heatmap_{uniprot_id}.png")
    plt.savefig(output_file) #save file 
    #plt.show() # Show plot

    plt.close() # Close the figure after saving to free up resources
