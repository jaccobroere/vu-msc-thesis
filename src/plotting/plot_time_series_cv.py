import os

import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use(["science"])
plt.rc("text", usetex=True)


def generate_matrices() -> tuple:
    """
    Generate two matrices of size p x p where each element is either an integer or 0.

    Parameters
    ----------
    p : int
        The size of the matrices.
    h : int
        The maximum distance between non-zero elements in the matrices.

    Returns
    -------
    tuple
        A tuple of two matrices of size p x p.

    """
    matrix = np.full((5, 15), 0.5, dtype=float)

    for i in range(5):
        matrix[i, i : i + 10] = 0.15
        matrix[i, i + 10] = 0.85

    return matrix


def plot_side_by_side_matrices() -> None:
    """
    Create a plot showing two matrices side by side with the non-zero elements numbered.

    Parameters
    ----------
    p : int
        The size of the matrices.
    h : int
        The maximum distance between non-zero elements in the matrices.

    Returns
    -------
    None

    """
    matrix = generate_matrices()

    fig, ax = plt.subplots(1, 1, figsize=(6, 2))
    ax.matshow(
        matrix,
        cmap="coolwarm",
        vmin=0,
        vmax=1,
    )

    for i in range(5):
        ax.axhline(i - 0.5, linestyle="-", color="k", linewidth=0.5)

    for i in range(15):
        ax.axvline(i - 0.5, linestyle="-", color="k", linewidth=0.5)

    # for j in range(9):
    #     ax.axhline(j - 0.5, linestyle="-", color="k", linewidth=0.5)
    #     ax.axvline(j - 0.5, linestyle="-", color="k", linewidth=0.5)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("Time")
    ax.set_ylabel("Folds")

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
    fig.savefig(filename)


if __name__ == "__main__":
    fig, ax = plot_side_by_side_matrices()
    save_figure(fig, os.path.join("out", "figures", "time_series_cv.eps"))
