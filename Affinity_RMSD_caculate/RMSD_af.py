#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates the RMSD between predicted AF2 structures (in PDBQT format)
and the corresponding real structures from the 'selected' folder,
based on residue ranges specified in 'benchmark_index.txt'.

Directory structure example:
pdbqt/
  af2_pdbqt/
    1e2k/
      1e2k.pdbqt
    1e2l/
      1e2l.pdbqt
    ...
selected/
  1e2k/
    1e2k_protein.pdb
  1e2l/
    1e2l_protein.pdb
  ...

The predicted structure is assumed to have CA atoms numbered from 0,
while the real structure uses actual residue numbering.
We use Bio.PDB.Superimposer to align the sets of CA atoms and compute RMSD.

Output lines will be in the format:
pdb_id <tab> RMSD_value
"""

import os
import re
from Bio import PDB

# ================ User Configuration ================
INDEX_FILE       = "benchmark_index.txt"    # Contains lines like: "1e2k  Chain A  Residues 55-60  length=6  DGPHGM"
AF2_PDBQT_DIR    = "pdbqt/af2_pdbqt"        # Directory containing subfolders for each pdb_id, each with a .pdbqt
REAL_STRUCT_DIR  = "selected"               # Directory containing real structure subfolders {pdbid}/{pdbid}_protein.pdb
OUTPUT_RMSD_FILE = "../Dataset_info/plt/with_af2/result_summary/af2_rmsd_summary.txt"  # Output file
# ====================================================

def parse_index_file(index_path):
    """
    Parse the benchmark index file, which contains lines like:
      1e2k  Chain A  Residues 55-60  length=6  DGPHGM
    Return a dictionary: { pdb_id: (chain_id, start_res, end_res) }.
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
            chain_id = parts[2]  # e.g. "A"
            match = residue_pat.search(line_s)
            if match:
                start_res = int(match.group(1))
                end_res   = int(match.group(2))
                index_dict[pdb_id] = (chain_id, start_res, end_res)
    return index_dict

def parse_pdbqt(pdbqt_path, length):
    """
    Parse the predicted AF2 PDBQT file to extract CA atom coordinates for the segment.
    The predicted structure is assumed to have residues numbered starting from 0.
    We only read up to 'length' CA atoms. Each CA line is expected to have x,y,z in columns 6,7,8 when split.
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
    Parse the real structure using Bio.PDB, returning CA coordinates
    from the specified chain for residues in [start_res, end_res].
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

def compute_rmsd(ref_coords, alt_coords):
    """
    Use Bio.PDB.Superimposer to align ref_coords and alt_coords, then return RMSD.
    Wrap the coordinate tuples with a dummy class providing get_coord().
    """
    from Bio.PDB import Superimposer

    if len(ref_coords) != len(alt_coords) or len(ref_coords) == 0:
        return None

    class DummyAtom:
        def __init__(self, coord):
            self.coord = coord
        def get_coord(self):
            return self.coord

    dummy_ref = [DummyAtom(c) for c in ref_coords]
    dummy_alt = [DummyAtom(c) for c in alt_coords]

    sup = Superimposer()
    sup.set_atoms(dummy_ref, dummy_alt)
    return sup.rms

def main():
    # Parse the index file to get chain and residue range
    index_dict = parse_index_file(INDEX_FILE)

    # Open output file
    with open(OUTPUT_RMSD_FILE, "w", encoding="utf-8") as outf:
        outf.write("# RMSD fig between AF2 .pdbqt and real structures\n\n")

        # For each pdb_id in index_dict, proceed
        for pdb_id, (chain_id, start_res, end_res) in index_dict.items():
            seg_length = end_res - start_res + 1
            if seg_length <= 0:
                outf.write(f"{pdb_id}\tN/A (invalid residue range)\n")
                continue

            # Predicted .pdbqt path: e.g. 'pdbqt/af2_pdbqt/1e2k/1e2k.pdbqt'
            pred_folder = os.path.join(AF2_PDBQT_DIR, pdb_id)
            pred_file   = os.path.join(pred_folder, f"{pdb_id}.pdbqt")
            if not os.path.isfile(pred_file):
                outf.write(f"{pdb_id}\tN/A (predicted file not found: {pred_file})\n")
                continue

            # Real PDB: e.g. 'selected/1e2k/1e2k_protein.pdb'
            real_folder = os.path.join(REAL_STRUCT_DIR, pdb_id)
            real_file   = os.path.join(real_folder, f"{pdb_id}_protein.pdb")
            if not os.path.isfile(real_file):
                outf.write(f"{pdb_id}\tN/A (real file not found: {real_file})\n")
                continue

            # Parse predicted coords (0..seg_length-1)
            pred_coords = parse_pdbqt(pred_file, seg_length)
            # Parse real coords (chain_id, start_res..end_res)
            real_coords = parse_real_pdb(real_file, chain_id, start_res, end_res)

            if len(pred_coords) != len(real_coords) or len(pred_coords) == 0:
                outf.write(f"{pdb_id}\tN/A (atom count mismatch: pred={len(pred_coords)}, real={len(real_coords)})\n")
                continue

            # Compute RMSD
            rmsd_val = compute_rmsd(real_coords, pred_coords)
            if rmsd_val is None:
                outf.write(f"{pdb_id}\tN/A (cannot compute RMSD)\n")
                continue

            # Write result in the format: 'pdb_id<TAB>RMSD_value'
            outf.write(f"{pdb_id}\t{rmsd_val:.3f}\n")
            print(f"{pdb_id} => RMSD={rmsd_val:.3f}")

    print(f"\nAll done! Results saved in {OUTPUT_RMSD_FILE}")

if __name__ == "__main__":
    main()
