# --*-- conding:utf-8 --*--
# @Time : 4/3/25 1:09 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : summary.py

"""
Generate a consolidated docking summary report by reading the final average
values from each protein's summary.txt file. Assumes each protein has its
own subdirectory under ROOT_DIR containing a summary.txt with a final average line.

- ROOT_DIR: folder containing one subdirectory per pdb_id
- SUMMARY_FILENAME: name of the per-protein summary file
- OUTPUT_SUMMARY_FILE: path to write the combined report
"""
import os
import re
from pathlib import Path

# Configuration
ROOT_DIR = "quantum_docking_result"  # Directory containing per-protein subfolders
SUMMARY_FILENAME = "summary.txt"     # Name of the summary file in each subfolder
OUTPUT_SUMMARY_FILE = os.path.join(
    "..", "results", "with_af3", "result_summary", "docking_summary_q.txt"
)  # File for the combined report

# Regex to extract final average values from a summary line
FINAL_PATTERN = re.compile(
    r"affinity\s*=\s*([\-\d\.]+).*rmsd_l\.b\.=\s*([\-\d\.]+).*rmsd_u\.b\.=\s*([\-\d\.]+)",
    re.IGNORECASE
)


def parse_final_average(summary_path):
    """
    Read a summary.txt and return a tuple (affinity, rmsd_lb, rmsd_ub).
    Expects a line starting with '[FINAL AVERAGE ACROSS' followed by a
    line matching the FINAL_PATTERN. Returns (None, None, None) if not found.
    """
    if not summary_path.is_file():
        return (None, None, None)

    found_final = False
    with summary_path.open('r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if found_final:
                match = FINAL_PATTERN.search(line)
                if match:
                    try:
                        aff = float(match.group(1))
                        lb  = float(match.group(2))
                        ub  = float(match.group(3))
                        return (aff, lb, ub)
                    except ValueError:
                        return (None, None, None)
            if line.startswith('[FINAL AVERAGE ACROSS'):
                found_final = True
    return (None, None, None)


def gather_summaries(root_dir):
    """
    Traverse each pdb_id subfolder under root_dir, parse its summary.txt,
    and collect final average values into a dict {pdb_id: (aff, lb, ub)}.
    """
    results = {}
    for entry in os.listdir(root_dir):
        protein_folder = Path(root_dir) / entry
        if not protein_folder.is_dir():
            continue
        summary_file = protein_folder / SUMMARY_FILENAME
        aff, lb, ub = parse_final_average(summary_file)
        if aff is not None:
            results[entry] = (aff, lb, ub)
    return results


def write_combined_report(data, output_file):
    """
    Write the combined docking summary to output_file in the format:
    # All proteins summary

    [pdb_id]
    affinity=..., rmsd_lb=..., rmsd_ub=...

    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w', encoding='utf-8') as f:
        f.write("# All proteins summary\n\n")
        for pdb_id in sorted(data.keys()):
            aff, lb, ub = data[pdb_id]
            f.write(f"[{pdb_id}]\n")
            f.write(f"affinity={aff:.3f}, rmsd_lb={lb:.3f}, rmsd_ub={ub:.3f}\n\n")


def main():
    summaries = gather_summaries(ROOT_DIR)
    write_combined_report(summaries, OUTPUT_SUMMARY_FILE)
    print(
        f"Combined summary for {len(summaries)} proteins written to '{OUTPUT_SUMMARY_FILE}'"
    )


if __name__ == '__main__':
    main()

