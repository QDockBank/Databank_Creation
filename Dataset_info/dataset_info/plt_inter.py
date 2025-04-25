#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads an interactions summary file where each line contains a
two-letter interaction (e.g. "KP") and a frequency count.
It builds a symmetric 20x20 matrix for all standard amino acids.
The heatmap is plotted using seaborn.heatmap with the "RdBu_r" colormap,
with linear normalization from 0 to 38. Each cell in the heatmap is annotated
with the corresponding data value.
All fonts are set to Arial.
The final heatmap is saved as an image file.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==================== User Configuration ====================
INPUT_FILE = "interactions_summary.txt"  # Input file (e.g., "KP    10")
OUTPUT_IMAGE = "aa_interaction.png"  # Output image file for the heatmap

# Predefined list of all 20 standard amino acid one-letter codes
ALL_AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
                   'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                   'W', 'Y']

FIGSIZE = (8, 8)  # Figure size, adjust as needed


# ==============================================================

def read_interactions(filepath):
    """
    Read the interactions summary file.
    Each valid line is expected to be in the format:
       INTERACTION<TAB>COUNT
    For example: "KP    10"
    Returns a list of tuples: [(aa1, aa2, count), ...]
    """
    interactions = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            interaction_str = parts[0]
            try:
                count = float(parts[1])
            except ValueError:
                continue
            if len(interaction_str) == 2:
                aa1, aa2 = interaction_str[0], interaction_str[1]
                interactions.append((aa1, aa2, count))
    return interactions


def build_interaction_matrix(interactions, all_aas):
    """
    Build a symmetric interaction matrix for all amino acids in all_aas.
    Rows and columns correspond to amino acids in the provided list.
    For each interaction (aa1, aa2, count), if both amino acids are in all_aas,
    add the count to the corresponding cell and ensure symmetry.

    Returns a numpy array of shape (n, n), where n = len(all_aas).
    """
    n = len(all_aas)
    matrix = np.zeros((n, n), dtype=float)
    aa_index = {aa: i for i, aa in enumerate(all_aas)}

    for aa1, aa2, count in interactions:
        if aa1 in aa_index and aa2 in aa_index:
            i = aa_index[aa1]
            j = aa_index[aa2]
            matrix[i, j] += count
            if i != j:
                matrix[j, i] += count
    return matrix


def plot_heatmap(all_aas, matrix, output_file):
    """
    Plot a heatmap of the interaction matrix using seaborn.
    Uses the "RdBu_r" colormap with linear normalization (vmin=0, vmax=38).
    Each cell is annotated with its numerical value (formatted to one decimal).
    The x-axis and y-axis display all standard amino acids.
    """
    # Set font to Arial
    plt.rcParams["font.family"] = "Arial"

    norm = plt.Normalize(vmin=0, vmax=17)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    # Plot heatmap with annotation enabled and specify the format (e.g., one decimal place)
    hm = sns.heatmap(matrix, cmap="RdBu_r", norm=norm,
                     annot=True, fmt=".0f", annot_kws={"fontsize": 10},
                     xticklabels=all_aas, yticklabels=all_aas, ax=ax,
                     linewidths=2, linecolor="white")

    # Rotate x-axis labels for clarity
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax.set_title("Amino Acid Interaction Frequency", fontsize=14)

    plt.tight_layout()
    plt.savefig(output_file, dpi=600)
    plt.show()
    print(f"Heatmap saved to {output_file}")


def main():
    interactions = read_interactions(INPUT_FILE)
    if not interactions:
        print("No interactions found in the input file.")
        return
    matrix = build_interaction_matrix(interactions, ALL_AMINO_ACIDS)
    plot_heatmap(ALL_AMINO_ACIDS, matrix, OUTPUT_IMAGE)


if __name__ == "__main__":
    main()
