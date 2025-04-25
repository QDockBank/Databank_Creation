# --*-- conding:utf-8 --*--
# @Time : 3/25/25 6:47 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Cifprepare.py

import os
import shutil
import subprocess
import numpy as np

from Bio.PDB import MMCIFParser, PDBIO, PDBParser


class Cifprepare:
    """
    A class to process protein structures from .cif files and convert them into .pdbqt files.

    Workflow:
      1) read_cif()            : Read the CIF file (using Biopython) and store structure in memory
      2) translate_to_origin() : (Optional) Translate the entire protein so that its center of mass is at the origin
      3) write_pdb()           : Write out the (optionally translated) structure to a PDB file
      4) prepare_pdbqt()       : Add hydrogens, compute partial (Gasteiger) charges, and convert to PDBQT via Open Babel
      5) run_pipeline()        : Run the above steps in sequence
    """

    def __init__(self, cif_file, output_folder=None):
        """
        Parameters:
        -----------
        cif_file : str
            Path to the .cif file containing the protein structure
        output_folder : str, optional
            The folder to save the final .pdbqt file. Default is the .cif .
        """
        self.cif_file = cif_file
        self.structure = None  # Biopython Structure object

        self.base_dir = os.path.dirname(os.path.abspath(cif_file))
        base_name = os.path.splitext(os.path.basename(cif_file))[0]


        self.output_pdb = os.path.join(self.base_dir, f"{base_name}.pdb")
        self.output_pdbqt = os.path.join(self.base_dir, f"{base_name}.pdbqt")

        if output_folder is None:
            self.docking_folder = self.base_dir
        else:
            self.docking_folder = output_folder
            os.makedirs(self.docking_folder, exist_ok=True)

        print(f"CIF base directory: {self.base_dir}")

    def _is_tool_available(self, tool_name):
        """Check if a specified tool (e.g., obabel) is available in the system."""
        return subprocess.call(f"type {tool_name}", shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

    def read_cif(self):
        """
        Read the CIF file using MMCIFParser (from Biopython) and store the result in self.structure.
        """
        if not os.path.exists(self.cif_file):
            raise FileNotFoundError(f"{self.cif_file} does not exist.")
        parser = MMCIFParser(QUIET=True)
        try:
            self.structure = parser.get_structure("cif_structure", self.cif_file)
            print(f"Successfully read CIF file: {self.cif_file}")
        except Exception as e:
            raise ValueError(f"Failed to parse CIF file {self.cif_file}: {e}")

    def calculate_center_of_mass(self):
        """
        Calculate the center of mass of the current protein structure stored in self.structure.
        Returns:
        --------
        center_of_mass : np.ndarray, shape (3,)
        """
        if self.structure is None:
            raise ValueError("Structure data is not loaded. Call read_cif() first.")

        coords = []
        for atom in self.structure.get_atoms():
            coords.append(atom.get_coord())
        if not coords:
            raise ValueError("No atoms found in the structure.")

        coords = np.array(coords)
        center = coords.mean(axis=0)
        print(f"Center of mass: X={center[0]:.3f}, Y={center[1]:.3f}, Z={center[2]:.3f}")
        return center

    def translate_to_origin(self):
        """
        Translate all atoms in the structure so that the center of mass is moved to the origin (0,0,0).
        """
        if self.structure is None:
            raise ValueError("Structure data is not loaded. Call read_cif() first.")

        center = self.calculate_center_of_mass()
        translation = -center

        for atom in self.structure.get_atoms():
            atom.set_coord(atom.get_coord() + translation)

        print("The structure has been translated so that its center of mass is at the origin.")

    def write_pdb(self, pdb_path=None):
        """
        Write the current structure to a PDB file.

        Parameters:
        -----------
        pdb_path : str, optional
            The path to the output PDB file. If None, use self.output_pdb.
        """
        if self.structure is None:
            raise ValueError("Structure data is not loaded. Call read_cif() first.")

        if pdb_path is None:
            pdb_path = self.output_pdb

        io = PDBIO()
        io.set_structure(self.structure)
        io.save(pdb_path)
        print(f"PDB file saved as: {pdb_path}")

    def prepare_pdbqt(self, translate=True, output_pdbqt=None):
        """
        Convert the current structure to a PDBQT file with hydrogens and Gasteiger charges using Open Babel.

        Parameters:
        -----------
        translate : bool
            Whether to translate the center of mass to the origin before conversion.
        output_pdbqt : str, optional
            The path to the output .pdbqt file. If None, use self.output_pdbqt.
        """
        if not self._is_tool_available("obabel"):
            raise EnvironmentError("Open Babel is not installed or not in PATH.")

        if translate:
            print("Translating to origin before generating PDBQT...")
            self.translate_to_origin()

        self.write_pdb(self.output_pdb)

        if output_pdbqt is None:
            output_pdbqt = self.output_pdbqt

        print("Converting to PDBQT using Open Babel...")
        command = [
            "obabel", self.output_pdb, "-O", output_pdbqt,
            "--partialcharge", "gasteiger", "-h", "-xr"
        ]
        subprocess.run(command, check=True)
        print(f"Protein PDBQT file has been generated: {output_pdbqt}")

        if os.path.abspath(self.docking_folder) != os.path.abspath(os.path.dirname(output_pdbqt)):
            final_pdbqt_path = os.path.join(self.docking_folder, os.path.basename(output_pdbqt))
            shutil.move(output_pdbqt, final_pdbqt_path)
            print(f"PDBQT file has been moved to: {final_pdbqt_path}")

    def run_pipeline(self, translate=True):
        """
        Run the pipeline from CIF -> PDB -> PDBQT (adding hydrogens and charges).

        Parameters:
        -----------
        translate : bool
            Whether to translate the center of mass to the origin before conversion.
        """
        print("==== 1) Reading CIF ====")
        self.read_cif()
        print("==== 2) Preparing PDBQT ====")
        self.prepare_pdbqt(translate=translate)
        print("==== Pipeline DONE ====")
