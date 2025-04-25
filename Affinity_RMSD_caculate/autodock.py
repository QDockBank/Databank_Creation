# --*-- conding:utf-8 --*--
# @Time : 3/24/25 2:44 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : autodock.py


import os
import random
from docking import AutoDockDocking


def main():

    protein_pdbqt_dir = "pdbqt/protein_pdbqt"
    ligand_pdbqt_dir = "pdbqt/ligand_pdbqt"
    docking_output_dir = "docking_results"
    os.makedirs(docking_output_dir, exist_ok=True)


    for protein_id in os.listdir(protein_pdbqt_dir):
        protein_subdir = os.path.join(protein_pdbqt_dir, protein_id)
        if not os.path.isdir(protein_subdir):
            continue

        ligand_subdir = os.path.join(ligand_pdbqt_dir, protein_id)
        ligand_file = os.path.join(ligand_subdir, f"{protein_id}_ligand.pdbqt")
        if not os.path.exists(ligand_file):
            print(f"[Warning] Ligand not found for {protein_id}, expected {ligand_file}, skipping.")
            continue


        seeds_for_runs = [random.randint(1, 1000000) for _ in range(20)]
        print(f"[Info] Protein {protein_id} => seeds: {seeds_for_runs}")


        fragment_files = [f for f in os.listdir(protein_subdir) if f.endswith(".pdbqt")]
        fragment_files.sort()  # e.g. 1bai_top_1.pdbqt, 1bai_top_2.pdbqt, ...

        for frag_file in fragment_files:
            frag_path = os.path.join(protein_subdir, frag_file)


            frag_name_no_ext = os.path.splitext(frag_file)[0]
            frag_output_dir = os.path.join(docking_output_dir, protein_id, frag_name_no_ext)
            os.makedirs(frag_output_dir, exist_ok=True)


            for run_idx in range(20):
                run_seed = seeds_for_runs[run_idx]

                run_output_subdir = os.path.join(frag_output_dir, f"run{run_idx + 1}")
                os.makedirs(run_output_subdir, exist_ok=True)

                docking_obj = AutoDockDocking(
                    receptor_pdbqt=frag_path,
                    ligand_pdbqt=ligand_file,
                    output_dir=run_output_subdir,
                    seed=run_seed
                )

                out_pdbqt, log_file = docking_obj.run_docking()

                seed_record_file = os.path.join(run_output_subdir, "seed_used.txt")
                with open(seed_record_file, "w") as sf:
                    sf.write(str(run_seed))

                print(f"[Done] {protein_id}/{frag_file} run {run_idx + 1} => {out_pdbqt}")

    print("\nAll docking processes completed!")

if __name__ == "__main__":
    main()

