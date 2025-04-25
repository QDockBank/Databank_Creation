# --*-- conding:utf-8 --*--
# @time:4/24/25 07:24
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File:create_metadata.py


"""
Split a single metadata TXT file into per‐protein JSON files.

Each JSON will have two top‐level keys:
  - "protein_information": { "pdb_id", "sequence", "sequence_length", "residues" }
  - "quantum_metadata": { "number_of_qubits", "circuit_depth",
                         "lowest_energy", "highest_energy",
                         "energy_range", "execution_time_s" }

Assumes:
  - The metadata TXT is tab‐separated, with a header row matching the example.
  - There is a parent dataset directory containing one subfolder per pdb_id.
"""

import os
import json
import pandas as pd

# Path to the metadata TXT file
METADATA_FILE = "dataset_info/benchmark_info.txt"

# Base directory under which each pdb_id has its own subfolder
BASE_DIR = "QDockBank"

def main():
    # Read the metadata table
    # df = pd.read_csv(METADATA_FILE, sep='\t', dtype=str)
    df = pd.read_csv(METADATA_FILE, sep='\t', skipinitialspace=True, dtype=str)

    # Clean up column names and map to JSON keys
    df = df.rename(columns={
        'pdb_id':             'pdb_id',
        'Residue_sequence':   'sequence',
        'Sequence_length':    'sequence_length',
        'Residues':           'residues',
        'Number_of_qubits':   'number_of_qubits',
        'Circuit_Depth':      'circuit_depth',
        'Lowest_energy':      'lowest_energy',
        'Highest_energy':     'highest_energy',
        'Energy_range':       'energy_range',
        'Execution_Time(s)':  'execution_time_s',
    })

    # Process each row and write out JSON
    for _, row in df.iterrows():
        pdb_id = row['pdb_id']
        folder = os.path.join(BASE_DIR, pdb_id)
        os.makedirs(folder, exist_ok=True)

        # Build the JSON structure
        record = {
            "protein_information": {
                "pdb_id":          pdb_id,
                "sequence":        row['sequence'],
                "sequence_length": int(row['sequence_length']),
                "chain": "A",
                "residues":        row['residues'],
            },
            "quantum_metadata": {
                "number_of_qubits": int(row['number_of_qubits']),
                "circuit_depth":    int(row['circuit_depth']),
                "lowest_energy":    float(row['lowest_energy']),
                "highest_energy":   float(row['highest_energy']),
                "energy_range":     float(row['energy_range']),
                "execution_time_s": float(row['execution_time_s']),
            }
        }

        # Write JSON to {pdb_id}/{pdb_id}_metadata.json
        out_path = os.path.join(folder, f"{pdb_id}_metadata.json")
        with open(out_path, 'w') as f:
            json.dump(record, f, indent=2)

        print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
