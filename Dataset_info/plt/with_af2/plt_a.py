# --*-- conding:utf-8 --*--
# @Time : 4/7/25 1:49 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : plt_a.py

import os
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.family'] = 'Arial'


def load_data(filename):
    data = []
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            pdb_id = parts[0]
            try:
                quantum_val = float(parts[1].split("=")[1])
                af2_val = float(parts[2].split("=")[1])
            except (IndexError, ValueError):
                continue
            better = parts[3].split("=")[1]
            data.append({
                "pdb_id": pdb_id,
                "quantum": quantum_val,
                "af2": af2_val,
                "better": better
            })
    df = pd.DataFrame(data)
    return df


def plot_scatter(df, save_path=None):
    # Split data by category ("better")
    df_quantum = df[df['better'] == "quantum"]
    df_af3 = df[df['better'] == "af2"]

    plt.figure(figsize=(6, 6))

    # Plot quantum points with circle marker and af3 points with triangle marker
    plt.scatter(df_quantum['quantum'], df_quantum['af2'], s=150, marker='o', color='tab:blue', edgecolors="k",
                linewidths=1, alpha=0.95, label="Better: Quantum")
    plt.scatter(df_af3['quantum'], df_af3['af2'], s=150, marker='*', color='tab:orange', edgecolors="k", linewidths=1,
                alpha=0.95, label="Better: AF2")

    x_min = df['quantum'].min() - 0.5
    x_max = df['quantum'].max() + 0.5
    y_min = df['af2'].min() - 0.5
    y_max = df['af2'].max() + 0.5
    lim_min = min(x_min, y_min)
    lim_max = max(x_max, y_max)
    plt.xlim(lim_min, lim_max)
    plt.ylim(lim_min, lim_max)

    plt.plot([lim_min, lim_max], [lim_min, lim_max], "k--", label="x = y")

    plt.xlabel("Quantum", fontsize=14, fontname="Arial")
    plt.ylabel("AF2", fontsize=14, fontname="Arial")
    plt.title("Compare Affinity: Quantum vs AF2", fontsize=16, fontname="Arial")
    plt.legend(fontsize=10, prop={'family': 'Arial'})

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format="png", dpi=600, bbox_inches="tight")
        print(f"Plot saved to {save_path}")
    plt.show()


def main():
    group = 'S'
    data_file = f"{group}/affinity_compare.txt"
    df = load_data(data_file)
    print(df)
    output_plot_folder = f"plots/{group}"
    os.makedirs(output_plot_folder, exist_ok=True)
    save_path = os.path.join(output_plot_folder, f"Affinity_af2_{group}.png")
    plot_scatter(df, save_path=save_path)


if __name__ == "__main__":
    main()