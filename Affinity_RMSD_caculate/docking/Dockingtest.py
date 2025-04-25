# --*-- conding:utf-8 --*--
# @Time : 11/10/24 3:40 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Dockingtest.py

import os
import subprocess


class AutoDockDocking:
    def __init__(self, receptor_pdbqt, ligand_pdbqt, output_dir, log_file_name=None, seed=2):
        """
        Parameters:
        -----------
        receptor_pdbqt : str
            Path to the receptor PDBQT file.
        ligand_pdbqt : str
            Path to the ligand PDBQT file (already converted).
        output_dir : str
            Directory to save docking fig.
        log_file_name : str, optional
            Log file name; if not provided, a dynamic name based on receptor and ligand will be used.
        seed : int, optional
            Random seed for docking.
        """
        self.receptor_pdbqt = receptor_pdbqt
        self.ligand_pdbqt = ligand_pdbqt
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.seed = seed

        # Dynamically generate a log file name if not provided.
        if log_file_name is None:
            receptor_base = os.path.splitext(os.path.basename(receptor_pdbqt))[0]
            ligand_base = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
            self.log_file_name = f"{receptor_base}_{ligand_base}_docking_log.txt"
        else:
            self.log_file_name = log_file_name

    def run_docking(self):
        """
        Run docking using AutoDock Vina with dynamic output file names.
        """
        # Use the provided ligand PDBQT file directly.
        # Calculate the center of mass of the ligand from its PDBQT file.
        center_x, center_y, center_z = self.calculate_center_of_mass(self.ligand_pdbqt)

        # Dynamically generate output file names based on receptor and ligand names.
        receptor_base = os.path.splitext(os.path.basename(self.receptor_pdbqt))[0]
        ligand_base = os.path.splitext(os.path.basename(self.ligand_pdbqt))[0]
        output_file = os.path.join(self.output_dir, f"docking_output_{receptor_base}_{ligand_base}.pdbqt")
        log_file = os.path.join(self.output_dir, self.log_file_name)

        print("Running AutoDock Vina docking...")
        self._run_vina(self.receptor_pdbqt, self.ligand_pdbqt, output_file, log_file,
                       center_x, center_y, center_z)

        print(f"Docking complete. Results saved in {output_file}")
        return output_file, log_file

    def calculate_center_of_mass(self, pdbqt_file):
        """
        Calculate the center of mass of a molecule from a PDBQT file.
        Assumes the file uses standard PDB format for ATOM records.
        """
        atom_coords = []
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        # PDB coordinates are typically in columns 31-38, 39-46, 47-54
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        atom_coords.append((x, y, z))
                    except ValueError:
                        continue
        if not atom_coords:
            raise ValueError(f"No atom coordinates found in file {pdbqt_file}.")
        xs, ys, zs = zip(*atom_coords)
        center_x = sum(xs) / len(xs)
        center_y = sum(ys) / len(ys)
        center_z = sum(zs) / len(zs)
        print(f"Calculated center of mass from {pdbqt_file}: X={center_x:.3f}, Y={center_y:.3f}, Z={center_z:.3f}")
        return center_x, center_y, center_z

    def _run_vina(self, receptor_pdbqt, ligand_pdbqt, output_file, log_file, center_x, center_y, center_z):
        """
        Run AutoDock Vina docking with the specified parameters.
        """
        # Set docking box size dynamically (可根据需要调整)
        size_x, size_y, size_z = 18, 18, 18
        command = [
            "vina",
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--out", output_file,
            "--center_x", str(center_x),
            "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--exhaustiveness", "16",
            "--seed", str(self.seed)
        ]
        try:
            with open(log_file, "w") as log:
                subprocess.run(command, stdout=log, stderr=log, check=True)
        except subprocess.CalledProcessError as e:
            print("Error: AutoDock Vina failed with the following command:")
            print(" ".join(command))
            print(f"Exit status: {e.returncode}")

    def _is_tool_available(self, tool_name):
        """Check if a tool is available on the system."""
        return subprocess.call(f"type {tool_name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

