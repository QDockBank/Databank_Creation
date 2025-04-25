# --*-- conding:utf-8 --*--
# @Time : 4/3/25 1:01 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : get_average_q.py


"""
Compute average docking statistics and write summary files for each protein.

This script scans a root folder containing one subfolder per protein (named by pdb_id).
Each protein folder contains multiple run subfolders (run1, run2, ...), each with a
ligand_docking_log.txt file. The script computes per-run averages of affinity and RMSD
(lower and upper bound), then computes the overall average across all runs,
writing a summary.txt in each protein folder.

Assumptions:
  - ROOT_DIR contains only pdb_id subfolders.
  - Each run folder is named 'runX' and contains exactly one file ending with '_docking_log.txt'.

Usage:
  python compute_docking_summary.py
"""

import os
import re
from pathlib import Path


ROOT_DIR = "quantum_docking_result"  # Root folder with one subfolder per pdb_id
SUMMARY_FILENAME = "summary.txt"     # Name of the summary file to write


# Pattern to match lines of the form:
#   mode_number   affinity_value   rmsd_lb   rmsd_ub
MODE_PATTERN = re.compile(r"^\s*(\d+)\s+([\-0-9\.]+)\s+([\-0-9\.]+)\s+([\-0-9\.]+)")


def parse_docking_log(log_file_path):
    """
    Parse a docking log file and compute the average affinity, lower RMSD, and upper RMSD.
    Return a tuple (affinity_mean, rmsd_lb_mean, rmsd_ub_mean) or (None, None, None) if no data.
    """
    total_aff, total_lb, total_ub = 0.0, 0.0, 0.0
    count = 0

    if not log_file_path.is_file():
        return (None, None, None)

    with log_file_path.open('r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            match = MODE_PATTERN.match(line)
            if not match:
                continue
            _, aff_str, lb_str, ub_str = match.groups()
            try:
                total_aff += float(aff_str)
                total_lb += float(lb_str)
                total_ub += float(ub_str)
                count += 1
            except ValueError:
                continue

    if count == 0:
        return (None, None, None)
    return (total_aff / count, total_lb / count, total_ub / count)


def compute_summary(protein_folder):
    """
    For a given protein folder, parse all run directories, compute per-run averages,
    then compute the overall average. Return a dict with 'per_run' list and 'final' tuple.
    """
    run_dirs = [d for d in protein_folder.iterdir() if d.is_dir() and d.name.startswith('run')]
    run_dirs.sort(key=lambda d: int(d.name.replace('run', '')) if d.name.replace('run', '').isdigit() else float('inf'))

    per_run_results = []  # list of (run_name, (aff, lb, ub))
    for run_dir in run_dirs:
        log_files = list(run_dir.glob('*_docking_log.txt'))
        if not log_files:
            continue
        # assume only one docking log per run
        aff, lb, ub = parse_docking_log(log_files[0])
        if aff is not None:
            per_run_results.append((run_dir.name, (aff, lb, ub)))

    if not per_run_results:
        return None

    # compute overall averages
    total_aff = sum(r[1][0] for r in per_run_results)
    total_lb  = sum(r[1][1] for r in per_run_results)
    total_ub  = sum(r[1][2] for r in per_run_results)
    n = len(per_run_results)

    final_aff = total_aff / n
    final_lb  = total_lb  / n
    final_ub  = total_ub  / n

    return {
        'per_run': per_run_results,
        'final': (final_aff, final_lb, final_ub)
    }


def write_summary(protein_folder, summary_data):
    """
    Write the summary.txt file in the protein_folder based on summary_data.
    """
    summary_path = protein_folder / SUMMARY_FILENAME
    with summary_path.open('w', encoding='utf-8') as f:
        f.write(f"# Docking Summary for {protein_folder.name}\n\n")
        for run_name, (aff, lb, ub) in summary_data['per_run']:
            f.write(f"{run_name} => affinity={aff:.3f}, rmsd_lb={lb:.3f}, rmsd_ub={ub:.3f}\n")
        aff_f, lb_f, ub_f = summary_data['final']
        f.write(f"\n# Final average over {len(summary_data['per_run'])} runs\n")
        f.write(f"affinity={aff_f:.3f}, rmsd_lb={lb_f:.3f}, rmsd_ub={ub_f:.3f}\n")


def main():
    root = Path(ROOT_DIR)
    if not root.is_dir():
        print(f"Error: ROOT_DIR does not exist: {ROOT_DIR}")
        return

    print(f"Scanning protein folders in '{ROOT_DIR}'...")
    for item in sorted(root.iterdir()):
        if not item.is_dir():
            continue
        summary = compute_summary(item)
        if summary is None:
            print(f"Warning: No valid run data for {item.name}")
            continue
        write_summary(item, summary)
        print(f"Summary written for {item.name}")

    print("All summaries generated.")


if __name__ == '__main__':
    main()

