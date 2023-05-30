import os

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["mathtext.fontset"] = "dejavuserif"


# %%
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
    matrix_1 = np.zeros((p, p), dtype=float)
    matrix_2 = np.zeros((p, p), dtype=float)

    for i in range(p):
        for j in range(p):
            if i == j:
                matrix_1[i, j] = 0.15
            elif abs(i - j) <= h:
                matrix_1[i, j] = 0.75 + abs(i - j) * 0.07
            else:
                matrix_1[i, j] = 0.15

        for j in range(p):
            if abs(i - j) <= h:
                matrix_2[i, j] = 0.75 + abs(i - j) * 0.07
            else:
                matrix_2[i, j] = 0.15

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
    print(matrix_1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(2 * max(5, p // 2), max(5, p // 2)))
    ax1.matshow(
        matrix_1,
        cmap="coolwarm",
        vmin=0,
        vmax=1,
    )
    ax2.matshow(
        matrix_2,
        cmap="coolwarm",
        vmin=0,
        vmax=1,
    )

    for i in range(p):
        for j in range(p):
            if matrix_1[i, j] > 0.5:
                if j < i:
                    direction_y = i - j
                    direction_x = j - i
                    ax1.quiver(
                        (i),
                        (j),
                        direction_x,
                        direction_y,
                        angles="xy",
                        scale_units="xy",
                        scale=1,
                        color="w",
                        zorder=5,
                        width=0.005,
                        headwidth=5,
                    )

            if matrix_2[i, j] > 0.5:
                if j < i:
                    direction_y = i - j
                    direction_x = j - i
                    ax2.quiver(
                        (i),
                        (j),
                        direction_x,
                        direction_y,
                        angles="xy",
                        scale_units="xy",
                        scale=1,
                        color="w",
                        zorder=5,
                        width=0.005,
                        headwidth=5,
                    )

        for ax in (ax1, ax2):
            ax.axhline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
            ax.axvline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
            ax.axhline(p - 0.5, linestyle="-", color="k", linewidth=0.5)
            ax.axvline(p - 0.5, linestyle="-", color="k", linewidth=0.5)

            ax.set_xticks([])
            ax.set_yticks([])

    plt.tight_layout(pad=2.5)
    plt.show()

    return fig, (ax1, ax2)


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
    p = 9
    h = (p - 1) // 4
    fig, (ax1, ax2) = plot_side_by_side_matrices(p, h)
    save_figure(fig, os.path.join("out", "figures", "symmetric_arrows_AB.eps"))