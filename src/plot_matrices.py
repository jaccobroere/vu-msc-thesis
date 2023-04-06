import matplotlib.pyplot as plt
import numpy as np


def plot_side_by_side_matrices(p, h):
    # Create the matrices
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

    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(2 * max(5, p // 2), max(5, p // 2)))
    ax1.matshow(matrix_1, cmap="binary", vmin=0, vmax=0)
    ax2.matshow(matrix_2, cmap="binary", vmin=0, vmax=0)
    # ax1.matshow(matrix_1, cmap='viridis')
    # ax2.matshow(matrix_2, cmap='viridis')

    # Add the element numbers
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
    # Display the plot
    plt.show()


# Example usage
plot_side_by_side_matrices(25, 6)
