import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots
from matplotlib.ticker import MultipleLocator, NullLocator

# Change to project directory
PROJ_DIR = os.environ["PROJ_DIR"]
os.chdir(PROJ_DIR)

# Set plotting style to academic
plt.style.use(["science", "ieee"])
plt.rc("text", usetex=True)


def load_data():
    df_mean = pd.read_csv(
        os.path.join("out", "simulation", "runtime", "runtime_mean.csv"),
        index_col=0,
    )
    df_std = pd.read_csv(
        os.path.join("out", "simulation", "runtime", "runtime_std.csv"),
        index_col=0,
    )

    matcher = np.vectorize(lambda x: re.search(r".*p(\d+).*", x).group(1))
    df_mean.index = matcher(df_mean.index.values)
    df_std.index = matcher(df_std.index.values)

    # Change column names
    df_mean.columns = [
        r"F-SPLASH($\lambda$)",
        r"SSF-SPLASH(0.5, $\lambda$)",
        r"SPLASH(0.5, $\lambda$)",
        r"GF-SPLASH(0.5, $\lambda$)",
        r"PVAR($\lambda$)",
    ]

    return df_mean, df_std


def plot_runtimes_linegraph():
    # Read the data
    df_mean, df_std = load_data()

    # Set up the figure
    fig, ax = plt.subplots(figsize=(5, 3))

    # Add a line for each model
    for col in df_mean.iloc[:, :-1]:
        ax.plot(df_mean.index, df_mean[col], marker="o", label=col, markersize=2.5)

    # Change y-axis to log scale
    ax.set_yscale("log")
    ax.grid(True)

    # Remove x-axis minor ticks
    ax.xaxis.set_minor_locator(NullLocator())

    # Add a legend
    ax.legend()

    # Set axis labels
    ax.set_xlabel(r"Problem dimension ($p$)")
    ax.set_ylabel(r"Runtime ($s$)")

    # Show plot
    plt.tight_layout()
    plt.show()

    return fig, ax


def save_figure(fig, filename):
    """
    Save a figure to a file.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    filename : str
        The name of the file to save the figure to.

    Returns
    -------
    None
    """
    fig.savefig(filename, dpi=1000)


if __name__ == "__main__":
    fig, ax = plot_runtimes_linegraph()
    save_figure(fig, os.path.join("out", "figures", "runtime_linegraph.eps"))
