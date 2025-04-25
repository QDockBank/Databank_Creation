# --*-- conding:utf-8 --*--
# @Time : 2/21/25 6:05 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : vqe_prediction.py

import os
import shutil
import time

from Protein_Folding import Peptide
from Protein_Folding.interactions.miyazawa_jernigan_interaction import MiyazawaJerniganInteraction
from Protein_Folding.penalty_parameters import PenaltyParameters
from Protein_Folding.protein_folding_problem import ProteinFoldingProblem
from qiskit_ibm_runtime import QiskitRuntimeService

from Qiskit_VQE import VQE
from Qiskit_VQE import StateCalculator

def read_config(file_path):
    """
    Read the config file (INSTANCE=xxx / TOKEN=xxx)
    """
    config = {}
    try:
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                if not line or "=" not in line:
                    continue
                key, value = line.split("=", 1)
                config[key.strip()] = value.strip()
    except Exception as e:
        print(f"Failed to read config file: {e}")
        return None
    return config


def parse_txt_file(txt_file_path):
    """
    Read each line from the specified TXT file and parse it into the following docking_workspace structure:
    [
      {
        'pdb_id': <str>,
        'pocket_file': <str>,
        'chain': <str>,
        'residues_range': <str>,
        'sequence': <str>   # single-letter sequence
      },
      ...
    ]
    """
    # Mapping from three-letter codes to one-letter codes
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    fragments = []
    with open(txt_file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Assumes tab-separated docking_workspace. Adjust if your docking_workspace uses spaces or other delimiters
            parts = line.split('\t')
            if len(parts) < 5:
                continue

            pdb_id = parts[0]
            pocket_file = parts[1]
            chain_info = parts[2]          # e.g. "Chain A"
            residue_info = parts[3]        # e.g. "Residues 192-203"
            seq_info_3letter = parts[4]    # e.g. "VAL-VAL-TYR-PRO-..."

            # Parse chain
            chain = chain_info.replace("Chain", "").strip()
            # Parse residue range
            residues_range = residue_info.replace("Residues", "").strip()

            # Convert three-letter codes to one-letter codes
            three_letter_codes = seq_info_3letter.split('-')
            single_letter_codes = []
            for code_3letter in three_letter_codes:
                code_3letter = code_3letter.strip().upper()
                if code_3letter in three_to_one:
                    single_letter_codes.append(three_to_one[code_3letter])
                else:
                    raise ValueError(f"Unknown residue code: {code_3letter}")

            seq_info_1letter = "".join(single_letter_codes)

            fragments.append({
                'pdb_id': pdb_id,
                'pocket_file': pocket_file,
                'chain': chain,
                'residues_range': residues_range,
                'sequence': seq_info_1letter
            })
    return fragments


def pick_unique_fragments(fragments, max_count=25):
    """
    Keep the first max_count fragments with different pdb_ids
    """
    selected = []
    used_pdb_ids = set()
    for frag in fragments:
        if frag['pdb_id'] not in used_pdb_ids:
            selected.append(frag)
            used_pdb_ids.add(frag['pdb_id'])
        if len(selected) >= max_count:
            break
    return selected


def run_vqe_for_fragment(frag, service, max_iter=150):
    """
    Use VQE to predict the protein structure corresponding to the given fragment
    """
    main_chain_sequence = frag['sequence']
    protein_id = frag['pdb_id']

    print(f"\n=== Processing protein {protein_id} ===")
    print(f"Residue sequence: {main_chain_sequence}")
    print(f"Sequence length: {len(main_chain_sequence)}")

    # 1. Construct the protein object (only main chain here, no side chains)
    side_chain_seq = ['' for _ in range(len(main_chain_sequence))]
    peptide = Peptide(main_chain_sequence, side_chain_seq)

    # 2. Define interaction and penalty parameters
    mj_interaction = MiyazawaJerniganInteraction()
    penalty_terms = PenaltyParameters(10, 10, 10)

    # 3. Build the protein folding problem and construct the Hamiltonian
    protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
    hamiltonian = protein_folding_problem.qubit_op()

    # Here, according to your needs, add 5 extra qubits (instead of the original +2 approach)
    qubits_num = hamiltonian.num_qubits + 5
    print(f"Number of qubits: {qubits_num}")

    # 4. Call VQE
    vqe_instance = VQE(
        service=service,
        hamiltonian=hamiltonian,
        min_qubit_num=qubits_num,
        maxiter=max_iter
    )

    # run_vqe() returns (energy_list, best_solution, ansatz, top_results)
    energy_list, final_solution, ansatz, top_results = vqe_instance.run_vqe()

    # 5. Post-processing and saving fig
    # Set output directory: one folder per pdb_id
    output_dir = f"result/{protein_id}"
    os.makedirs(output_dir, exist_ok=True)

    # Save the entire iteration energy list
    energy_list_file = os.path.join(output_dir, f"energy_list_{protein_id}.txt")
    with open(energy_list_file, 'w') as file:
        for item in energy_list:
            file.write(str(item) + '\n')
    print(f"Energy list saved to: {energy_list_file}")

    # Compute the probability distribution of the final_solution and interpret it into 3D coordinates
    state_calculator = StateCalculator(service, qubits_num, ansatz)
    final_prob_dist = state_calculator.get_probability_distribution(final_solution)
    protein_result = protein_folding_problem.interpret(final_prob_dist)

    # Save the final structure as .xyz
    protein_result.save_xyz_file(name=protein_id, path=output_dir)
    print(f"Protein structure saved at: {output_dir}/{protein_id}.xyz")

    print(f"Finished processing: {protein_id}\n")


def main():

    txt_file_path = "Data/predicted_data.txt"  # The TXT file containing fragments to be predicted
    config_path = "config.txt"  # IBM Quantum config file
    max_fragments = 22                   # Maximum number of fragments
    max_iter = 200                       # Maximum VQE iterations

    config = read_config(config_path)
    if not config or "INSTANCE" not in config or "TOKEN" not in config:
        print("Could not read INSTANCE or TOKEN from config. Please check your config file.")
        return

    # Initialize the quantum service

    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance=config["INSTANCE"],
        token=config["TOKEN"]
    )

    # Read all fragments from TXT
    all_fragments = parse_txt_file(txt_file_path)

    # Pick up to max_fragments unique fragments
    selected_fragments = pick_unique_fragments(all_fragments, max_fragments)

    # Run quantum prediction for each fragment
    log_file_path = "execution_time_log.txt"
    with open(log_file_path, 'w') as log_file:
        log_file.write("Protein_ID\tSequence\tExecution_Time(s)\n")

        for idx, fragment in enumerate(selected_fragments, start=1):
            protein_id = fragment['pdb_id']
            seq = fragment['sequence']
            print(f"\n--- Fragment {idx}/{len(selected_fragments)}: PDB={protein_id} ---")

            start_time = time.time()
            run_vqe_for_fragment(fragment, service, max_iter=max_iter)
            end_time = time.time()

            elapsed = end_time - start_time
            # Record to log
            log_file.write(f"{protein_id}\t{seq}\t{elapsed:.2f}\n")

    print("\nAll processing is complete. Log saved to:", log_file_path)


if __name__ == "__main__":
    main()


