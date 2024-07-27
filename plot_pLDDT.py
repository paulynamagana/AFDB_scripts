import matplotlib.pyplot as plt

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
