#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script: plot raw data bar charts (quantum, af2, af3) without normalization or scaling.
Only annotate the quantum bars with their values; af2 and af3 bars are not annotated.
Colors and bar widths are user-configurable.

Example file structure:
  - with_af2/all/affinity_compare.txt: line format e.g. "1abc  quantum=-5.12  af2=-4.98  better=af2"
  - with_af3/all/affinity_compare.txt: similar, for quantum vs af3
  - with_af2/all/rmsd_compare.txt, with_af3/all/rmsd_compare.txt: similar, for RMSD
Adjust file paths, input methods, etc. as needed.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# ==================== User Configurable Section ====================

# --- File paths ---
AFF_FILE_AF2 = "with_af2/all/affinity_compare.txt"
AFF_FILE_AF3 = "with_af3/all/affinity_compare.txt"

RMSD_FILE_AF2 = "with_af2/all/rmsd_compare.txt"
RMSD_FILE_AF3 = "with_af3/all/rmsd_compare.txt"

# --- Output image files ---
OUTPUT_AFF_PLOT = "fig/affinity_plot.png"
OUTPUT_RMSD_PLOT = "fig/rmsd_plot.png"

# --- Color and bar width settings ---
AFF_COLOR_QUANTUM = "cyan"
AFF_COLOR_AF2 = "lightgray"
AFF_COLOR_AF3 = "gray"

RMSD_COLOR_QUANTUM = "greenyellow"
RMSD_COLOR_AF2 = "lightgray"
RMSD_COLOR_AF3 = "gray"

BAR_WIDTH_QUANTUM = 0.2
BAR_WIDTH_AF2 = 0.2
BAR_WIDTH_AF3 = 0.2

# --- Figure sizes ---
AFF_FIG_SIZE = (13, 3)
RMSD_FIG_SIZE = (13, 3)

# Set global font to Arial
plt.rcParams["font.family"] = "Arial"

# ==================== Function Implementations ====================

def read_compare_file(filepath):
    """
    Read a comparison file and return a dict { pdb_id: (quantum_val, method_val) }.
    Line examples:
      1abc  quantum=-5.123  af2=-4.567  better=quantum
      or
      1abc  quantum=-3.12   af3=-3.20   better=af3
    """
    data = {}
    if not os.path.isfile(filepath):
        print(f"[Warning] File not found: {filepath}")
        return data

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            pdb_id = parts[0]
            try:
                quantum_val = float(parts[1].split("=")[1])
                method_val = float(parts[2].split("=")[1])
                data[pdb_id] = (quantum_val, method_val)
            except Exception as e:
                print(f"Error parsing line '{line}': {e}")
    return data


def merge_comparison_data(file_af2, file_af3):
    """
    Take two files (quantum vs af2 and quantum vs af3) and return a dict { pdb_id: (quantum, af2, af3) }.
    """
    data_af2 = read_compare_file(file_af2)
    data_af3 = read_compare_file(file_af3)

    merged = {}
    all_ids = set(data_af2.keys()) | set(data_af3.keys())
    for pdb_id in all_ids:
        q_val = None
        af2_val = None
        af3_val = None
        if pdb_id in data_af2:
            q_val, af2_val = data_af2[pdb_id]
        if pdb_id in data_af3:
            q2_val, temp_af3 = data_af3[pdb_id]
            if q_val is None:
                q_val = q2_val
            af3_val = temp_af3
        if q_val is not None:
            merged[pdb_id] = (q_val, af2_val, af3_val)
    return merged


def plot_bar_chart(metric_data, output_file, fig_size, ylabel,
                   color_quantum, color_af2, color_af3,
                   w_q, w_a2, w_a3):
    """
    Plot grouped bar charts (quantum, af2, af3) with actual values.
    Only annotate the quantum bars; af2/af3 bars are not annotated.
    """
    pdb_ids = sorted(metric_data.keys())
    if not pdb_ids:
        print(f"No data to plot for {ylabel}.")
        return

    val_q = []
    val_af2 = []
    val_af3 = []
    valid_ids = []

    for pid in pdb_ids:
        q, a2, a3 = metric_data[pid]
        if q is None and a2 is None and a3 is None:
            continue
        valid_ids.append(pid)
        val_q.append(q if q is not None else 0)
        val_af2.append(a2 if a2 is not None else 0)
        val_af3.append(a3 if a3 is not None else 0)

    n = len(valid_ids)
    if n == 0:
        print(f"No valid entries for {ylabel}.")
        return

    x = np.arange(n)
    offset_q = -0.2
    offset_a2 = 0
    offset_a3 = 0.2

    fig, ax = plt.subplots(figsize=fig_size)

    # Plot quantum bars
    bars_q = ax.bar(x + offset_q, val_q, width=w_q,
                    color=color_quantum, edgecolor="black",
                    linewidth=0.5, label="quantum")
    # Plot af2 bars
    bars_af2 = ax.bar(x + offset_a2, val_af2, width=w_a2,
                      color=color_af2, edgecolor="black",
                      linewidth=0.5, label="af2")
    # Plot af3 bars
    bars_af3 = ax.bar(x + offset_a3, val_af3, width=w_a3,
                      color=color_af3, edgecolor="black",
                      linewidth=0.5, label="af3")

    # # Only annotate quantum bars
    # for i, bar in enumerate(bars_q):
    #     height = bar.get_height()
    #     ax.text(
    #         bar.get_x() + bar.get_width() / 2,
    #         height,
    #         f"{val_q[i]:.2f}",
    #         ha="center",
    #         va="bottom",
    #         rotation=90,
    #         fontsize=8
    #     )

    ax.set_xticks(x)
    ax.set_xticklabels(valid_ids, rotation=45, ha="right", fontsize=9)

    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f"{ylabel} Comparison (No Normalization)", fontsize=13)
    ax.legend(fontsize=10)
    ax.axhline(0, color="black", linewidth=1)

    plt.tight_layout()
    plt.savefig(output_file, dpi=600)
    plt.show()
    print(f"Plot saved => {output_file}")


def main():
    # 1. Merge raw data
    aff_data = merge_comparison_data(AFF_FILE_AF2, AFF_FILE_AF3)
    rmsd_data = merge_comparison_data(RMSD_FILE_AF2, RMSD_FILE_AF3)

    print(f"[Info] Affinity data size = {len(aff_data)}")
    print(f"[Info] RMSD data size = {len(rmsd_data)}")

    # 2. Plot Affinity chart (quantum, af2, af3)
    plot_bar_chart(
        aff_data,
        OUTPUT_AFF_PLOT,
        AFF_FIG_SIZE,
        ylabel="Affinity Value",
        color_quantum=AFF_COLOR_QUANTUM,
        color_af2=AFF_COLOR_AF2,
        color_af3=AFF_COLOR_AF3,
        w_q=BAR_WIDTH_QUANTUM,
        w_a2=BAR_WIDTH_AF2,
        w_a3=BAR_WIDTH_AF3
    )

    # 3. Plot RMSD chart (quantum, af2, af3)
    plot_bar_chart(
        rmsd_data,
        OUTPUT_RMSD_PLOT,
        RMSD_FIG_SIZE,
        ylabel="RMSD Value",
        color_quantum=RMSD_COLOR_QUANTUM,
        color_af2=RMSD_COLOR_AF2,
        color_af3=RMSD_COLOR_AF3,
        w_q=BAR_WIDTH_QUANTUM,
        w_a2=BAR_WIDTH_AF2,
        w_a3=BAR_WIDTH_AF3
    )


if __name__ == "__main__":
    main()
