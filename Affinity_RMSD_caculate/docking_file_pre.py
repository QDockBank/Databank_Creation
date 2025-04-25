# --*-- conding:utf-8 --*--
# @Time : 3/20/25 8:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : docking_file_pre.py

"""
Process protein and ligand XYZ/MOL2 files into PDBQT format.
Assumes each protein folder under grouped_result_dir/chain_x/ contains exactly one XYZ file named "{pdb_id}.xyz".
Ligand files are expected in PDBbind/{pdb_id}/{pdb_id}_ligand.mol2".
"""
import os
from docking import Fileprepare

# Configurable directories
grouped_result_dir = "grouped_prediction_result"
selected_dir = "PDBbind"
output_protein_dir = "pdbqt/protein_pdbqt"
output_ligand_dir = "pdbqt/ligand_pdbqt"

# Ensure output directories exist
os.makedirs(output_protein_dir, exist_ok=True)
os.makedirs(output_ligand_dir, exist_ok=True)

if __name__ == '__main__':
    for chain_dir in os.listdir(grouped_result_dir):
        chain_path = os.path.join(grouped_result_dir, chain_dir)
        if not os.path.isdir(chain_path):
            continue

        # Iterate each protein folder (named by pdb_id)
        for pdb_id in os.listdir(chain_path):
            protein_path = os.path.join(chain_path, pdb_id)
            if not os.path.isdir(protein_path):
                continue

            # Expect a single XYZ: {pdb_id}.xyz
            xyz_filename = f"{pdb_id}.xyz"
            xyz_path = os.path.join(protein_path, xyz_filename)
            if not os.path.isfile(xyz_path):
                print(f"Warning: Expected {xyz_filename} in {protein_path}, skipping.")
                continue

            # Process protein
            protein_output_dir = os.path.join(output_protein_dir, pdb_id)
            os.makedirs(protein_output_dir, exist_ok=True)
            print(f"\n=== Processing protein '{pdb_id}' ===")
            docking_obj = Fileprepare(
                xyz_file=xyz_path,
                docking_folder=protein_output_dir
            )
            docking_obj.run_pipeline()

            # Process ligand if exists
            ligand_mol2 = os.path.join(selected_dir, pdb_id, f"{pdb_id}_ligand.mol2")
            if os.path.isfile(ligand_mol2):
                ligand_output_dir = os.path.join(output_ligand_dir, pdb_id)
                os.makedirs(ligand_output_dir, exist_ok=True)
                print(f"=== Processing ligand for '{pdb_id}' ===")
                docking_obj.prepare_ligand_pdbqt(
                    ligand_mol2_path=ligand_mol2,
                    output_folder=ligand_output_dir
                )
            else:
                print(f"Warning: ligand file not found for {pdb_id} at {ligand_mol2}")

    print("\nAll processing done!")


