#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates the RMSD between the predicted top segment (from pdbqt files)
and the corresponding segment of the real PDB structure.
It reads a top-list file (with lines like "1e2k  top_1") and a benchmark index file,
then for each protein:
  - It retrieves the predicted pdbqt file specified by the top selection,
  - Extracts CA atom coordinates from the predicted pdbqt file for the segment (assuming predicted residue numbering starts at 0),
  - Extracts CA atom coordinates from the real PDB (using the specified residue range),
  - Uses Bio.PDB.Superimposer to superimpose the two sets of atoms and compute the RMSD,
  - Finally writes the result in a two-column format: "pdb_id<TAB>RMSD".
"""

import os
import re
from Bio import PDB

# ====================== User Configuration ======================
TOP_LIST_FILE = "create_benchmark/top_selected.txt"       # File containing lines like "1e2k  top_1"
INDEX_FILE = "benchmark_index.txt"                         # File containing chain and residue range, e.g. "1e2k  Chain A  Residues 55-60 ..."
PREDICTED_DIR = "pdbqt/protein_pdbqt"                      # Folder containing subfolders with predicted top .pdbqt files
REAL_STRUCT_DIR = "selected"                               # Folder containing real structure subfolders (each containing {pdb_id}_protein.pdb)
OUTPUT_RMSD_FILE = "../Dataset_info/plt/with_af3/result_summary/q_rmsd_summary.txt"  # Output file for RMSD fig (final format: pdb_id<TAB>RMSD)
# ===============================================================

def parse_top_list(filename):
    """
    Parse the file that contains lines in the format:
      pdb_id    top_name
    Return a dictionary mapping pdb_id to top_name.
    """
    top_dict = {}
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            if not line_s or line_s.startswith("#"):
                continue
            parts = line_s.split()
            if len(parts) < 2:
                continue
            pdb_id = parts[0]
            top_name = parts[1]
            top_dict[pdb_id] = top_name
    return top_dict

def parse_index_file(index_path):
    """
    Parse the benchmark index file that contains lines like:
      1e2k    Chain A    Residues 55-60    length=6    DGPHGM
    Return a dictionary mapping pdb_id to a tuple (chain_id, start_res, end_res).
    Skip lines that start with '[' or 'length_' or are blank.
    """
    index_dict = {}
    residue_pat = re.compile(r'Residues\s+(\d+)-(\d+)', re.IGNORECASE)
    with open(index_path, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            if not line_s:
                continue
            if line_s.startswith("[") or line_s.startswith("length_"):
                continue
            parts = line_s.split()
            if len(parts) < 5:
                continue
            pdb_id = parts[0]
            chain_id = parts[2]  # e.g., "A"
            match_r = residue_pat.search(line_s)
            if match_r:
                start_res = int(match_r.group(1))
                end_res = int(match_r.group(2))
                index_dict[pdb_id] = (chain_id, start_res, end_res)
    return index_dict

def parse_pdbqt(pdbqt_path, length):
    """
    Parse the predicted pdbqt file and extract CA atom coordinates for the segment.
    The predicted file is assumed to have residue numbering starting from 0.
    Only lines with ' CA ' are processed and the function stops when 'length' CA atoms are read.
    Returns a list of (x, y, z) tuples.
    """
    coords = []
    current_ca_count = 0
    with open(pdbqt_path, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            if not line_s.startswith("ATOM"):
                continue
            if " CA " not in line_s:
                continue
            parts = line_s.split()
            try:
                x = float(parts[6])
                y = float(parts[7])
                z = float(parts[8])
            except (IndexError, ValueError):
                continue
            coords.append((x, y, z))
            current_ca_count += 1
            if current_ca_count >= length:
                break
    return coords

def parse_real_pdb(pdb_path, chain_id, start_res, end_res):
    """
    Parse the real PDB structure using Bio.PDB and extract CA atom coordinates
    from the specified chain within the residue range [start_res, end_res].
    Returns a list of (x, y, z) tuples.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("real_struct", pdb_path)
    coords = []
    for model in structure:
        for chain in model:
            if chain.get_id() == chain_id:
                for residue in chain:
                    res_id = residue.get_id()[1]
                    if start_res <= res_id <= end_res:
                        if "CA" in residue:
                            atom = residue["CA"]
                            coords.append(tuple(atom.coord))
    return coords

def compute_rmsd(ref_atoms, alt_atoms):
    """
    Use Bio.PDB.Superimposer to superimpose ref_atoms and alt_atoms and return the RMSD.
    Both ref_atoms and alt_atoms are lists of coordinate tuples.
    We wrap each tuple in a DummyAtom that provides a get_coord() method.
    """
    from Bio.PDB import Superimposer

    if len(ref_atoms) != len(alt_atoms) or len(ref_atoms) == 0:
        return None

    # Define a simple DummyAtom class
    class DummyAtom:
        def __init__(self, coord):
            self.coord = coord
        def get_coord(self):
            return self.coord

    dummy_ref = [DummyAtom(coord) for coord in ref_atoms]
    dummy_alt = [DummyAtom(coord) for coord in alt_atoms]

    sup = Superimposer()
    sup.set_atoms(dummy_ref, dummy_alt)
    return sup.rms

def main():
    # Parse the top list file
    top_dict = parse_top_list(TOP_LIST_FILE)
    # Parse the benchmark index file for chain and residue range information
    index_dict = parse_index_file(INDEX_FILE)

    with open(OUTPUT_RMSD_FILE, "w", encoding="utf-8") as outf:
        # Write header comment (optional)
        outf.write("# RMSD fig (predicted vs real) in format: pdb_id <tab> RMSD_value\n\n")
        for pdb_id, top_name in top_dict.items():
            if pdb_id not in index_dict:
                outf.write(f"{pdb_id}\tN/A  (not in benchmark index)\n")
                continue
            chain_id, start_res, end_res = index_dict[pdb_id]
            seg_length = end_res - start_res + 1
            if seg_length <= 0:
                outf.write(f"{pdb_id}\tN/A  (invalid residue range {start_res}-{end_res})\n")
                continue

            # Predicted pdbqt file path (assumed naming: {pdb_id}_{top_name}.pdbqt in subfolder {pdb_id})
            pred_folder = os.path.join(PREDICTED_DIR, pdb_id)
            pred_filename = f"{pdb_id}_{top_name}.pdbqt"
            pred_path = os.path.join(pred_folder, pred_filename)
            if not os.path.isfile(pred_path):
                outf.write(f"{pdb_id}\tN/A  (predicted file not found: {pred_path})\n")
                continue

            # Real PDB file path (assumed naming: {pdb_id}_protein.pdb in subfolder {pdb_id})
            real_folder = os.path.join(REAL_STRUCT_DIR, pdb_id)
            real_pdb = os.path.join(real_folder, f"{pdb_id}_protein.pdb")
            if not os.path.isfile(real_pdb):
                outf.write(f"{pdb_id}\tN/A  (real PDB not found: {real_pdb})\n")
                continue

            # Parse predicted CA coordinates (using relative numbering: [0 .. seg_length-1])
            pred_coords = parse_pdbqt(pred_path, seg_length)
            # Parse real CA coordinates (using actual numbering: [start_res .. end_res])
            real_coords = parse_real_pdb(real_pdb, chain_id, start_res, end_res)

            if len(pred_coords) != len(real_coords) or len(pred_coords) == 0:
                outf.write(f"{pdb_id}\tN/A  (atom count mismatch: pred={len(pred_coords)}, real={len(real_coords)})\n")
                continue

            # Compute RMSD using Superimposer on the wrapped DummyAtom objects
            rmsd_val = compute_rmsd(pred_coords, real_coords)
            if rmsd_val is None:
                outf.write(f"{pdb_id}\tN/A  (cannot compute RMSD)\n")
                continue

            outf.write(f"{pdb_id}\t{rmsd_val:.3f}\n")
            print(f"{pdb_id} => RMSD={rmsd_val:.3f}")

    print(f"\nAll done! RMSD fig saved in {OUTPUT_RMSD_FILE}")

if __name__ == "__main__":
    main()
