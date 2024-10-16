import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_plddt_legend(dpi=100):
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
    plt.title('Average AM and pLDDT scores per position')
    plt.legend()
    plt.grid(axis='y', linestyle='--')  # Optional grid

        # Save the plot to a file
    output_directory = "data_output"
    
    output_file = os.path.join(output_directory, f"graph_plDDT-AM-score_{uniprot_id}.png")
    plt.savefig(output_file) #save file 
    
    # Show the plot
    plt.show()
    plt.close()

