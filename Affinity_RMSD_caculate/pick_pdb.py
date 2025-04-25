# --*-- conding:utf-8 --*--
# @time:4/9/25 10:48 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : pick_pdb.py

"""
This script creates a benchmark dataset by copying per-protein PDB files
into a standardized folder structure. It uses an index file to select
which pdb_ids to include and searches recursively in the source folder
for each "{pdb_id}_full_model_translated.pdb" file.

Assumptions:
  - index_group.txt lists pdb_ids under [Group S/M/L] sections.
  - Each pdb_id desired has exactly one matching PDB file in the source tree.
  - GROUPED_FOLDER and BENCHMARK_FOLDER are sibling directories.

Usage:
  python pick_pdb.py
"""

import os
import re
import shutil


INDEX_FILE = "index_group.txt"       # File listing pdb_ids under [Group S/M/L]
SOURCE_FOLDER = "grouped_prediction_result"  # Root folder to search for PDBs
BENCHMARK_FOLDER = "QDockbank" # Destination for cleaned PDBs


def parse_group_index(index_file):
    """
    Parse the index file, extracting all pdb_ids under [Group S/M/L] sections.
    """
    pdb_set = set()
    current_group = None
    with open(index_file, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            g = re.match(r'^\[Group\s+([SML])\]', s)
            if g:
                current_group = g.group(1)
                continue
            if current_group in {"S","M","L"}:
                if s.startswith("length_"):
                    continue
                pdb_set.add(s.split()[0])
    return pdb_set


def find_pdb_file(pdb_id, root_folder):
    """
    Recursively search for "{pdb_id}_full_model_translated.pdb" under root_folder.
    Returns the full path if found, else None.
    """
    target = f"{pdb_id}_full_model_translated.pdb"
    for dirpath, _, filenames in os.walk(root_folder):
        if target in filenames:
            return os.path.join(dirpath, target)
    return None


def main():
    pdb_set = parse_group_index(INDEX_FILE)
    print(f"Found {len(pdb_set)} pdb_ids in '{INDEX_FILE}'")

    os.makedirs(BENCHMARK_FOLDER, exist_ok=True)
    copied = 0

    for pdb_id in sorted(pdb_set):
        path = find_pdb_file(pdb_id, SOURCE_FOLDER)
        if not path:
            print(f"Warning: PDB file for {pdb_id} not found")
            continue

        dest_dir = os.path.join(BENCHMARK_FOLDER, pdb_id)
        os.makedirs(dest_dir, exist_ok=True)
        dest_file = os.path.join(dest_dir, f"{pdb_id}.pdb")

        shutil.copy2(path, dest_file)
        print(f"Copied {path} → {dest_file}")
        copied += 1

    print(f"\nDone: {copied} files copied into '{BENCHMARK_FOLDER}'")

if __name__ == "__main__":
    main()