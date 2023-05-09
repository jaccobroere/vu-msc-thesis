import matplotlib.pyplot as plt
import numpy as np


def generate_matrices(p: int, h: int) -> tuple:
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
    matrix_1 = np.zeros((p, p), dtype=int)
    matrix_2 = np.zeros((p, p), dtype=int)

    counter = 1
    for i in range(p):
        for j in range(p):
            if i == j:
                matrix_1[i, j] = 0
            elif abs(i - j) <= h:
                matrix_1[i, j] = counter
                counter += 1
            else:
                matrix_1[i, j] = 0

        for j in range(p):
            if abs(i - j) <= h:
                matrix_2[i, j] = counter
                counter += 1
            else:
                matrix_2[i, j] = 0

    return matrix_1, matrix_2


def plot_side_by_side_matrices(p: int, h: int) -> None:
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
    matrix_1, matrix_2 = generate_matrices(p, h)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(2 * max(5, p // 2), max(5, p // 2)))
    ax1.matshow(matrix_1, cmap="binary", vmin=0, vmax=0)
    ax2.matshow(matrix_2, cmap="binary", vmin=0, vmax=0)

    for i in range(p):
        for j in range(p):
            if matrix_1[i, j] != 0:
                ax1.text(j, i, f"{matrix_1[i, j]}", ha="center", va="center", color="k")
            if matrix_2[i, j] != 0:
                ax2.text(j, i, f"{matrix_2[i, j]}", ha="center", va="center", color="k")

        for ax in (ax1, ax2):
            ax.axhline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
            ax.axvline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
            ax.axhline(p - 0.5, linestyle="-", color="k", linewidth=0.5)
            ax.axvline(p - 0.5, linestyle="-", color="k", linewidth=0.5)

            ax.set_xticks([])
            ax.set_yticks([])

    plt.show()

    return fig


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


p = 9
h = p // 4
plot_side_by_side_matrices(p, h)
