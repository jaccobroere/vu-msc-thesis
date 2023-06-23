import os

import matplotlib.pyplot as plt
import numpy as np
import scienceplots
from matplotlib.colors import ListedColormap
from scipy.io import mmread

plt.style.use(["science"])
plt.rc("text", usetex=True)
plt.rcParams["text.latex.preamble"] = r"\usepackage{bm} \usepackage{amsmath}"


# %%
def my_cmap():
    colors = ["gainsboro", "red", "blue"]
    cmap = ListedColormap(colors)
    return cmap


# %%
def load_Dtilde(path: str) -> tuple:
    """
    Load the Dtilde matrix from a file and convert it to a numpy array.
    """
    Dtilde = mmread(path)

    return Dtilde.toarray()


# %%
def generate_matrices_color(matrix) -> tuple:
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
    matrix = matrix.copy()

    # 0 stays at 0
    # 1 stays at 1
    # Replace -1 values with third color
    matrix[matrix == -1] = 2
    # Replace alpha values with fourth color
    # Replace (1 - alpha) values with fifth color

    # matrix[matrix < 0.0] = 2
    # matrix[matrix > 0.0] = 1

    return matrix


# %%
def plot_penalty_matrices_FSPLASH_SSFSPLASH(dir_path: str) -> None:
    # Load the penalty matrices
    Dtilde_FSPLASH = generate_matrices_color(
        load_Dtilde(os.path.join(dir_path, "Dtilde.mtx"))
    )
    Dtilde_SSFSPLASH = generate_matrices_color(
        load_Dtilde(os.path.join(dir_path, "Dtilde_SSF.mtx"))
    )

    # Load colormap
    cmap = my_cmap()

    # Initialize the figure
    fig, axs = plt.subplots(1, 2, figsize=(6, 3))

    # Plot the matrices
    axs[0].matshow(
        Dtilde_FSPLASH,
        cmap=cmap,
        vmin=0,
        vmax=2,
    )
    axs[1].matshow(
        Dtilde_SSFSPLASH,
        cmap=cmap,
        vmin=0,
        vmax=2,
    )

    # for ax in axs.flatten():
    #     ax.axhline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
    #     ax.axvline(i - 0.5, linestyle="-", color="k", linewidth=0.5)
    #     ax.axhline(p - 0.5, linestyle="-", color="k", linewidth=0.5)
    #     ax.axvline(p - 0.5, linestyle="-", color="k", linewidth=0.5)

    #     ax.set_xticks([])
    #     ax.set_yticks([])

    # Set titles of the matrices
    axs[0].set_title(f"$\\tilde{{\\bm{{D}}}}^{{\\text{{F-SPLASH}}}}$", fontsize=12)
    axs[1].set_title(
        f"$\\tilde{{\\bm{{D}}}}^{{\\text{{SSF-SPLASH}}}}_{{\\alpha}}$", fontsize=12
    )

    max_y = max(Dtilde_FSPLASH.shape[0], Dtilde_SSFSPLASH.shape[0])

    for ax in axs.flatten():
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim(max_y - 0.5, -0.5)

    # Set the spacing between the subplots
    plt.tight_layout(pad=1.5)
    plt.show()

    return fig, axs


# %%
def plot_penalty_matrices_GFSPLASH(dir_path: str) -> None:
    # Load the penalty matrices
    D_alpha = generate_matrices_color(
        load_Dtilde(os.path.join(dir_path, "D_alpha_GFSPLASH.mtx"))
    )
    D_alpha_sym = generate_matrices_color(
        load_Dtilde(os.path.join(dir_path, "D_alpha_sym_GFSPLASH.mtx"))
    )

    res = []
    for row in D_alpha_sym:
        if not any((D_alpha == row).all(1)):
            res.append(row)

    sym_rows = np.array(res)

    Dtilde_alpha_sym = np.vstack((D_alpha, sym_rows))

    # Load colormap
    cmap = my_cmap()

    # Initialize the figure
    fig, axs = plt.subplots(
        1, 2, figsize=(6, 3 * D_alpha_sym.shape[0] / D_alpha.shape[1])
    )

    # Plot the matrices
    axs[0].matshow(
        D_alpha,
        cmap=cmap,
        vmin=0,
        vmax=2,
    )
    axs[1].matshow(
        Dtilde_alpha_sym,
        cmap=cmap,
        vmin=0,
        vmax=2,
    )

    # Set titles of the matrices
    axs[0].set_title(
        f"$\\tilde{{\\bm{{D}}}}^{{\\text{{GF-SPLASH}}}}_{{\\alpha, 0}}$", fontsize=12
    )
    axs[1].set_title(
        f"$\\tilde{{\\bm{{D}}}}^{{\\text{{GF-SPLASH}}}}_{{\\alpha, 1}}$", fontsize=12
    )

    for ax in axs.flatten():
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")

    # Set the spacing between the subplots
    plt.tight_layout(pad=1.5)
    plt.show()

    return fig, axs


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
    dir_path = os.path.join("out", "figures", "penalty_matrices")
    # Plot the penalty matrices for FSPLASH and SSFSPLASH
    fig, axs = plot_penalty_matrices_FSPLASH_SSFSPLASH(dir_path=dir_path)
    save_figure(
        fig, os.path.join("out", "figures", "penalty_matrices_FSPLASH_SSFSPLASH.eps")
    )
    # Plot the penalty matrices for GFSPLASH
    fig, axs = plot_penalty_matrices_GFSPLASH(dir_path=dir_path)
    save_figure(fig, os.path.join("out", "figures", "penalty_matrices_GFSPLASH.eps"))
