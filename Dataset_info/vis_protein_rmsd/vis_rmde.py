# --*-- conding:utf-8 --*--
# @time:4/13/25 6:04 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : benchmark.py

"""
This script manually extracts a specified fragment from a reference (pdbbind) pdb file and a predicted pdb file.
It then aligns the predicted fragment to the reference fragment using CA atoms and saves the aligned predicted fragment
as a new pdb file.
"""

import os
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer

# ================== User-specified parameters ==================
# Manual specification of target pdb id and chain
PDB_ID = "2qbs"  # Example: "1e1x"
CHAIN_ID = "A"  # Chain identifier

# Residue range in the reference (pdbbind) pdb file (standard structure)
REF_START_RES = 214  # e.g., start at residue 47
REF_END_RES = 224  # e.g., end at residue 59

# For the predicted structure, assume fragment is renumbered starting at 1.
# The predicted fragment length should be the same as reference fragment length.
PRED_START = 1
PRED_END = REF_END_RES - REF_START_RES + 1

# Directories and file paths
DATASET_DIR = "dataset"  # Directory containing standard pdb files, e.g. "1e1x_protein.pdb"
PREDICTED_FILE = '2qbs.pdb'  # Predicted pdb file for the given pdb id
OUTPUT_FILE = f"aligned_fragment_{PDB_ID}.pdb"  # Output file name for aligned predicted fragment


# ===============================================================

class ResidueSelect(Select):
    """
    A Select class for PDBIO to save only residues in a given residue number range from a specified chain.
    """

    def __init__(self, chain_id, start_res, end_res):
        self.chain_id = chain_id
        self.start_res = start_res
        self.end_res = end_res

    def accept_residue(self, residue):
        # The residue id is a tuple, second element is the residue number
        resnum = residue.get_id()[1]
        if self.start_res <= resnum <= self.end_res:
            return 1
        return 0


def extract_ca_atoms(structure, chain_id, start_res, end_res):
    """
    Extract CA atoms from the first model of the structure in the specified chain and residue range.
    Returns a list of CA atoms.
    """
    ca_atoms = []
    model = structure[0]
    chain = model[chain_id]
    for residue in chain:
        resnum = residue.get_id()[1]
        if start_res <= resnum <= end_res:
            if "CA" in residue:
                ca_atoms.append(residue["CA"])
    return ca_atoms


def extract_fragment_structure(structure, chain_id, start_res, end_res):
    """
    Extract a fragment substructure (only residues in the specified range from the chain)
    from the given structure.
    Returns a new structure containing only the pdbbind residues.
    """
    io = PDBIO()
    io.set_structure(structure)
    temp_filename = "temp_fragment.pdb"
    io.save(temp_filename, select=ResidueSelect(chain_id, start_res, end_res))

    parser = PDBParser(QUIET=True)
    fragment_struct = parser.get_structure("fragment", temp_filename)
    os.remove(temp_filename)  # Remove temporary file
    return fragment_struct


def main():
    parser = PDBParser(QUIET=True)

    # Construct paths for the standard and predicted pdb files
    standard_pdb = os.path.join(DATASET_DIR, f"{PDB_ID}_protein.pdb")
    if not os.path.isfile(standard_pdb):
        print(f"Reference pdb file not found: {standard_pdb}")
        return
    if not os.path.isfile(PREDICTED_FILE):
        print(f"Predicted pdb file not found: {PREDICTED_FILE}")
        return

    # Parse the structures
    ref_structure = parser.get_structure(f"ref_{PDB_ID}", standard_pdb)
    pred_structure = parser.get_structure(f"pred_{PDB_ID}", PREDICTED_FILE)

    # Extract CA atoms for the reference fragment from the standard pdb
    ref_ca_atoms = extract_ca_atoms(ref_structure, CHAIN_ID, REF_START_RES, REF_END_RES)
    if not ref_ca_atoms:
        print("No CA atoms found in reference fragment.")
        return

    # For the predicted structure, assume the corresponding fragment is renumbered starting at 1.
    pred_ca_atoms = extract_ca_atoms(pred_structure, CHAIN_ID, PRED_START, PRED_END)
    if not pred_ca_atoms:
        print("No CA atoms found in predicted fragment.")
        return

    if len(ref_ca_atoms) != len(pred_ca_atoms):
        print(f"Mismatch in number of CA atoms: reference {len(ref_ca_atoms)} vs predicted {len(pred_ca_atoms)}")
        return

    # Use Superimposer to align the predicted fragment to the reference fragment
    sup = Superimposer()
    sup.set_atoms(ref_ca_atoms, pred_ca_atoms)
    print(f"Alignment RMSD: {sup.rms:.3f}")
    sup.apply(pred_structure.get_atoms())  # Apply transformation to the entire predicted structure if desired

    # Extract the aligned predicted fragment as a new structure
    aligned_fragment = extract_fragment_structure(pred_structure, CHAIN_ID, PRED_START, PRED_END)

    # Save the aligned predicted fragment to OUTPUT_FILE
    io = PDBIO()
    io.set_structure(aligned_fragment)
    io.save(OUTPUT_FILE)
    print(f"Aligned fragment saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
