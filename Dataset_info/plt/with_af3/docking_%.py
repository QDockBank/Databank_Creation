#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script compares the affinity values between quantum and af2 docking fig.
The input files (for quantum and af2) now have the new format, where each line is:
    pdb_id<TAB>value
For example:
    1e2k    -3.138
The script reads the group index file (with [Group S/M/L]) and the two docking summary files,
then produces an overall comparison file and group-specific files.
The output format for each line is:
    pdb_id    quantum=-3.138    af2=-2.752    better=quantum
(Where lower value is better.)
"""

import os
import re

# ====================== User Configuration ======================
INDEX_FILE = "group_index.txt"  # Group index file, containing [Group S/M/L] and protein IDs.
QUANTUM_FILE = "result_summary/docking_summary_q.txt"  # Quantum docking fig file (new format).
AF2_FILE = "result_summary/docking_summary_af3.txt"  # AF2 docking fig file (new format).

# Overall output file
AFF_ALL_FILE = "all/affinity_compare.txt"

# Group-specific output files for S, M, L groups
S_AFF_FILE = "S/affinity_compare.txt"
M_AFF_FILE = "M/affinity_compare.txt"
L_AFF_FILE = "L/affinity_compare.txt"


# ==============================================================

def parse_group_index(index_file):
    """
    Parse the group index file (with [Group S/M/L]) and return a dictionary:
    {
      "S": set([pdb_id, ...]),
      "M": set([...]),
      "L": set([...])
    }
    """
    groups_map = {"S": set(), "M": set(), "L": set()}
    current_group = None
    with open(index_file, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            if not line_s:
                continue
            # Match "[Group X]"
            grp_match = re.match(r'^\[Group\s+([SML])\]', line_s)
            if grp_match:
                current_group = grp_match.group(1)
                continue
            # Skip lines starting with "length_"
            if line_s.startswith("length_"):
                continue
            # Otherwise treat the first token as pdb_id
            if current_group in groups_map:
                pdb_id = line_s.split()[0]
                groups_map[current_group].add(pdb_id)
    return groups_map


def parse_docking_file(file_path):
    """
    Parse the new-format docking result file.
    Each line is expected to be in the format:
        pdb_id<TAB>value
    Returns a dictionary: { pdb_id: value }
    """
    data_map = {}
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            # Skip empty lines and comments
            if not line_s or line_s.startswith("#"):
                continue
            parts = line_s.split()
            if len(parts) < 2:
                continue
            pdb_id = parts[0]
            try:
                value = float(parts[1])
                data_map[pdb_id] = value
            except ValueError:
                continue
    return data_map


def compare_affinity(quantum_data, af2_data):
    """
    Compare affinity values between quantum and af2.
    Lower value is better.
    Returns a list of tuples:
      (pdb_id, quantum_value, af2_value, better, line_str)
    where line_str is formatted as:
      "pdb_id    quantum=-3.138    af2=-2.752    better=quantum"
    Only pdb_ids existing in both datasets are compared.
    """
    results = []
    common_ids = set(quantum_data.keys()) & set(af2_data.keys())
    for pdb_id in sorted(common_ids):
        q_val = quantum_data[pdb_id]
        a_val = af2_data[pdb_id]
        if q_val < a_val:
            better = "quantum"
        elif q_val > a_val:
            better = "af3"
        else:
            better = "tie"
        line_str = f"{pdb_id}\tquantum={q_val:.3f}\taf3={a_val:.3f}\tbetter={better}"
        results.append((pdb_id, q_val, a_val, better, line_str))
    return results


def write_compare_all(results, out_path):
    """
    Write the overall affinity comparison fig to out_path.
    """
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# Compare affinity: quantum vs af3 (lower is better)\n\n")
        for pdb_id, q_val, a_val, better, line_str in results:
            f.write(line_str + "\n")
        total = len(results)
        if total == 0:
            f.write("\n# No data.\n\n")
            return
        quantum_better = sum(1 for _, _, _, b, _ in results if b == "quantum")
        af2_better = sum(1 for _, _, _, b, _ in results if b == "af3")
        tie_count = total - quantum_better - af2_better
        qb_pct = quantum_better / total * 100
        ab_pct = af2_better / total * 100
        tie_pct = tie_count / total * 100
        f.write(f"\n# total={total}\n")
        f.write(f"# quantum better: {quantum_better} ({qb_pct:.1f}%)\n")
        f.write(f"# af2 better: {af2_better} ({ab_pct:.1f}%)\n")
        if tie_count > 0:
            f.write(f"# tie: {tie_count} ({tie_pct:.1f}%)\n")
        f.write("\n")


def distribute_and_stats(results, group_map, s_file, m_file, l_file):
    """
    Distribute the comparison fig into group-specific files based on group_map.
    """
    S_list, M_list, L_list = [], [], []
    for pdb_id, q_val, a_val, better, line_str in results:
        if pdb_id in group_map["S"]:
            S_list.append((pdb_id, q_val, a_val, better, line_str))
        elif pdb_id in group_map["M"]:
            M_list.append((pdb_id, q_val, a_val, better, line_str))
        elif pdb_id in group_map["L"]:
            L_list.append((pdb_id, q_val, a_val, better, line_str))

    def write_list_and_stats(rows, f_out):
        quantum_better = sum(1 for _, _, _, b, _ in rows if b == "quantum")
        af2_better = sum(1 for _, _, _, b, _ in rows if b == "af3")
        tie_count = len(rows) - quantum_better - af2_better
        for _, _, _, _, line_str in rows:
            f_out.write(line_str + "\n")
        total = len(rows)
        if total == 0:
            f_out.write("\n# No data.\n\n")
            return
        qb_pct = quantum_better / total * 100
        ab_pct = af2_better / total * 100
        tie_pct = tie_count / total * 100
        f_out.write(f"\n# total={total}\n")
        f_out.write(f"# quantum better: {quantum_better} ({qb_pct:.1f}%)\n")
        f_out.write(f"# af3 better: {af2_better} ({ab_pct:.1f}%)\n")
        if tie_count > 0:
            f_out.write(f"# tie: {tie_count} ({tie_pct:.1f}%)\n")
        f_out.write("\n")

    write_list_and_stats(S_list, s_file)
    write_list_and_stats(M_list, m_file)
    write_list_and_stats(L_list, l_file)


def main():
    # Parse group index file to obtain groups map
    group_map = parse_group_index(INDEX_FILE)
    print("Parsed group index:", group_map)

    # Load docking fig for quantum and af2 methods (new format)
    quantum_data = parse_docking_file(QUANTUM_FILE)
    af2_data = parse_docking_file(AF2_FILE)
    print("Quantum entries:", len(quantum_data))
    print("AF3 entries:", len(af2_data))

    # Compare affinity values between quantum and af2
    affinity_results = compare_affinity(quantum_data, af2_data)
    print("Affinity comparison entries:", len(affinity_results))

    # Write overall affinity comparison file
    write_compare_all(affinity_results, AFF_ALL_FILE)

    # Write group-specific comparison files
    with open(S_AFF_FILE, "w", encoding="utf-8") as fs, \
            open(M_AFF_FILE, "w", encoding="utf-8") as fm, \
            open(L_AFF_FILE, "w", encoding="utf-8") as fl:
        fs.write("# S group - affinity compare\n\n")
        fm.write("# M group - affinity compare\n\n")
        fl.write("# L group - affinity compare\n\n")
        distribute_and_stats(affinity_results, group_map, fs, fm, fl)

    print("Affinity comparison completed!")
    print("Overall and group-specific affinity comparison files generated.")

if __name__ == "__main__":
    main()
