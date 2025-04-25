# --*-- conding:utf-8 --*--
# @time:4/13/25 4:10 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : benchmark.py

"""
This script reads a benchmark txt file that contains columns including:
  pdb_id, Residue_sequence, Sequence_length, Number_of_qubits, Top1_best_energy,
  Cost_max, Cost_range, Execution_Time(s), Circuit_Depth.
It then extracts the 'Residue_sequence' for each protein fragment and analyzes
the inter-residue interactions. For each fragment, it considers all pairs (i, j) with i < j.
If a pair (first_letter + second_letter) or its reverse has been seen before in that fragment,
it will not be counted again. Finally, the script aggregates the counts of these interactions
across all fragments and writes the results to an output txt file.
"""

import os
import sys
from collections import Counter


INPUT_FILE = "benchmark_info.txt"  # Input file containing benchmark data (tab-delimited)
OUTPUT_FILE = "interactions_summary.txt"  # Output file to save aggregated interactions


# ============================================================

def parse_benchmark_file(filepath):
    """
    Parse the benchmark file and extract a list of residue sequences.
    Assumes the file is tab-delimited with a header. The column "Residue_sequence"
    is expected to be the second column.

    Returns:
        A list of residue sequences (strings).
    """
    sequences = []
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if not lines:
        print("Input file is empty.")
        return sequences

    # Assume first line is header; find index of "Residue_sequence" (case-insensitive)
    header = lines[0].strip().split("\t")
    try:
        seq_index = next(i for i, col in enumerate(header) if col.lower() == "residue_sequence")
    except StopIteration:
        # If not found, assume it is the second column (index 1)
        seq_index = 1

    # Process data lines
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) <= seq_index:
            continue
        sequence = cols[seq_index].strip()
        if sequence:
            sequences.append(sequence)
    return sequences


def get_interactions_from_sequence(seq):
    """
    For a given residue sequence (string), compute all unique unordered interactions.
    For every pair (i, j) with i < j, form the candidate as seq[i] + seq[j].
    If either candidate or its reverse is already recorded in this fragment, skip it.

    Returns:
        A set of interaction strings.
    """
    interactions = set()
    n = len(seq)
    for i in range(n - 1):
        for j in range(i + 1, n):
            candidate = seq[i] + seq[j]
            reverse_candidate = seq[j] + seq[i]
            # If either candidate or its reverse is already present, skip adding candidate
            if candidate in interactions or reverse_candidate in interactions:
                continue
            interactions.add(candidate)
    return interactions


def aggregate_interactions(sequences):
    """
    Aggregate interactions from a list of residue sequences.

    Returns:
        A Counter object mapping each interaction (string) to its count across all sequences.
    """
    interaction_counter = Counter()
    for seq in sequences:
        inter_set = get_interactions_from_sequence(seq)
        interaction_counter.update(inter_set)
    return interaction_counter


def write_interactions(counter, output_file):
    """
    Write the aggregated interaction counts to an output file.
    Each line will be in the format: interaction<TAB>count
    """
    with open(output_file, "w", encoding="utf-8") as f:
        for interaction, count in counter.items():
            f.write(f"{interaction}\t{count}\n")
    print(f"Interactions summary saved to '{output_file}'")


def main():
    sequences = parse_benchmark_file(INPUT_FILE)
    if not sequences:
        print("No sequences found, exiting.")
        sys.exit(1)
    print(f"Found {len(sequences)} sequences in the benchmark file.")

    interaction_counter = aggregate_interactions(sequences)
    print(f"Found {len(interaction_counter)} unique interactions across all sequences.")

    write_interactions(interaction_counter, OUTPUT_FILE)


if __name__ == "__main__":
    main()
