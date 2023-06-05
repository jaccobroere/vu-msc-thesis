import os

import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use(["science"])
plt.rc("text", usetex=True)
plt.rcParams["text.latex.preamble"] = r"\usepackage{bm}"


# %%
def generate_matrix_1_color(m: int):
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
    matrix_1 = np.full(shape=(m, m), fill_value=20, dtype=float)

    return matrix_1


def generate_matrix_1_text(m: int):
    matrix_1_text = np.zeros(shape=(m, m), dtype=int)
    counter = 1
    for i in range(m):
        for j in range(m):
            matrix_1_text[i, j] = counter
            counter += 1

    return matrix_1_text


def generate_matrix_2_color(m: int):
    p = m**2
    matrix_2_color = np.full(shape=(p, p), fill_value=20, dtype=float)

    for i in range(1, p + 1):
        for j in range(1, p + 1):
            if i == j:
                continue

            # Vertical neighbours
            if abs(i - j) == m:
                matrix_2_color[i - 1, j - 1] = 5

            # Horizontal neighbours
            not_side_edge = not (
                (i % m == 0)
                and ((j - 1) % m == 0)
                or ((j % m == 0) and ((i - 1) % m == 0))
            )
            if abs(i - j) == 1 and not_side_edge:
                matrix_2_color[i - 1, j - 1] = 1

            # Diagonal neighbours
            not_top_edge = (i > m) and (j > m)
            not_bottom_edge = (i <= (p - m)) and (j <= (p - m))
            diag_value = 10

            # Bottom right diagonal neighbour
            if not_side_edge and not_bottom_edge:
                if abs(i - j) == (m + 1):
                    matrix_2_color[i - 1, j - 1] = diag_value

            # Top right diagonal neighbour
            if not_side_edge and not_top_edge:
                if abs(i - j) == (m - 1):
                    matrix_2_color[i - 1, j - 1] = diag_value

            # Bottom left diagonal neighbour
            if not_bottom_edge and not_side_edge:
                if abs(i - j) == (m - 1):
                    matrix_2_color[i - 1, j - 1] = diag_value

            # Top left diagonal neighbour
            if not_top_edge and not_side_edge:
                if abs(i - j) == (m + 1):
                    matrix_2_color[i - 1, j - 1] = diag_value

    return matrix_2_color


def plot_side_by_side_matrices(m: int, p: int, h: int) -> None:
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
    matrix_1_color = generate_matrix_1_color(m)
    matrix_1_text = generate_matrix_1_text(m)
    matrix_2_color = generate_matrix_2_color(m)
    # matrix_color_1, matrix_color_2 = np.where(matrix_1 == 0, 0.15, 0.85), np.where(
    #     matrix_2 == 0, 0.15, 0.85
    # )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 3))
    ax1.matshow(
        matrix_1_color,
        cmap="tab20c",
        vmin=1,
        vmax=20,
    )
    ax2.matshow(
        matrix_2_color,
        cmap="tab20c",
        vmin=1,
        vmax=20,
    )
    for i in range(m):
        for j in range(m):
            ax1.text(
                j,
                i,
                f"$y_{{{matrix_1_text[i, j]}}}$",
                ha="center",
                va="center",
                color="k",
                fontsize=10,
            )

        ax1.axhline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
        ax1.axvline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
        ax1.axhline(m - 0.5, linestyle="-", color="k", linewidth=0.5)
        ax1.axvline(m - 0.5, linestyle="-", color="k", linewidth=0.5)

        ax1.set_xticks([])
        ax1.set_yticks([])

    for i in range(p):
        for j in range(p):
            pass

        ax2.axhline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
        ax2.axvline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
        ax2.axhline(p - 0.5, linestyle="-", color="k", linewidth=0.5)
        ax2.axvline(p - 0.5, linestyle="-", color="k", linewidth=0.5)

        ax2.set_xticks([])
        ax2.set_yticks([])

    # ax1.set_title(r"$\bm{A}$", fontsize=12)
    # ax2.set_title(r"$\bm{B}$", fontsize=12)

    plt.tight_layout(pad=1.5)
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
    m = 5
    p = m**2
    h = (p - 1) // 4
    fig, (ax1, ax2) = plot_side_by_side_matrices(m, p, h)
    save_figure(fig, os.path.join("out", "figures", "fusion_intuition.eps"))
