#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads an affinity comparison file and splits its content into
two separate files: one for 'quantum' values and one for 'af2' values.
Each output file will contain lines in the format:
    pdb_id<TAB>value
"""

import os

# User configuration
INPUT_FILE = "all/compare_affinity.txt"  # Input file containing lines in the specified format
OUTPUT_QUANTUM_FILE = "result_summary/docking_summary_q.txt"  # Output file for quantum values
OUTPUT_AF2_FILE = "result_summary/docking_summary_af2.txt"  # Output file for af2 values


def split_affinity_file(input_file, output_quantum, output_af2):
    """
    Read the input file and split its content into two files.

    Each valid data line is expected to have the format:
      pdb_id<TAB>quantum=-3.138<TAB>af2=-2.752<TAB>better=...

    This function extracts the numeric values after "quantum=" and "af2="
    and writes two files in the format:
      pdb_id<TAB>value
    """
    quantum_lines = []
    af2_lines = []

    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            # Skip empty or comment lines
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            pdb_id = parts[0]
            try:
                quantum_val = parts[1].split("=")[1]
                af2_val = parts[2].split("=")[1]
            except IndexError:
                continue
            quantum_lines.append(f"{pdb_id}\t{quantum_val}")
            af2_lines.append(f"{pdb_id}\t{af2_val}")

    with open(output_quantum, "w", encoding="utf-8") as fq:
        fq.write("\n".join(quantum_lines) + "\n")
    with open(output_af2, "w", encoding="utf-8") as fa:
        fa.write("\n".join(af2_lines) + "\n")


def main():
    split_affinity_file(INPUT_FILE, OUTPUT_QUANTUM_FILE, OUTPUT_AF2_FILE)
    print(f"Quantum affinity data saved to {OUTPUT_QUANTUM_FILE}")
    print(f"AF2 affinity data saved to {OUTPUT_AF2_FILE}")


if __name__ == "__main__":
    main()
