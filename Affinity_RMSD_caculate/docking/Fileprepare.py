# --*-- coding:utf-8 --*--
# @Time : 3/21/25 12:57 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Fileprepare.py


import os
import numpy as np
import subprocess
import shutil
from modeller import environ, log
from modeller.automodel import automodel, assess
from modeller import alignment, model
from Bio.PDB import PDBParser, PDBIO, MMCIFParser


class Fileprepare:
    """
    Use MODELLER to build a full-atom PDB from an XYZ file containing only Cα atoms,
    add hydrogens, balance charges, and convert to a PDBQT file.
    The final PDBQT file will be saved to the specified folder.

    Workflow:
      1) read_xyz()           : Read the XYZ file (containing only Cα coordinates and one-letter residue codes)
      2) adjust_scale()       : Scale coordinates so that adjacent Cα distance is 3.8 Å
      3) write_ca_pdb()       : Write the Cα-only PDB file (as a MODELLER template)
      4) prepare_alignment()  : Generate an alignment file for modeling
      5) generate_full_model(): Use automodel to build the full-atom model (complete PDB)
      6) prepare_pdbqt()      : Add hydrogens, balance charges, and convert to PDBQT file
      7) prepare_ligand_pdbqt(): Process specified MOL2 file, translate molecule center to origin,
                                 and convert to PDBQT file (for docking)
      8) run_pipeline()       : Run the above steps in sequence (for proteins)
    """

    def __init__(self, xyz_file, docking_folder=None):
        """
        Parameters:
        -----------
        xyz_file : str
            Path to the XYZ file containing only Cα coordinates and one-letter residue codes
        docking_folder : str, optional
            Specify the folder to save the final PDBQT file, default is the xyz file's directory
        """
        self.xyz_file = xyz_file
        self.sequence = []
        self.coordinates = []
        self.scaled_coordinates = []

        self.base_dir = os.path.dirname(os.path.abspath(self.xyz_file))
        print("XYZ base directory:", self.base_dir)

        base_name = os.path.splitext(os.path.basename(self.xyz_file))[0]  # e.g. "1fkf_top_3"

        self.ca_pdb_filename = os.path.join(self.base_dir, f"{base_name}_ca.pdb")
        self.alignment_filename = os.path.join(self.base_dir, f"{base_name}.ali")
        self.full_model_filename = os.path.join(self.base_dir, f"{base_name}_full_model.pdb")
        self.pdbqt_filename = os.path.join(self.base_dir, f"{base_name}.pdbqt")

        self.template_code = os.path.splitext(os.path.basename(self.ca_pdb_filename))[0]

        if docking_folder is None:
            self.docking_folder = self.base_dir
        else:
            self.docking_folder = docking_folder
            os.makedirs(self.docking_folder, exist_ok=True)

        self.chain_id = 'A'

    def read_xyz(self):
        """Read the XYZ file, extracting one-letter residue codes and Cα coordinates. Skip the first two lines."""
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()[2:]  # Skip the first two lines (typically atom count and comment)
            for line in lines:
                parts = line.strip().split()
                if len(parts) == 4:
                    res_name = parts[0]
                    x, y, z = map(float, parts[1:])
                    self.sequence.append(res_name)
                    self.coordinates.append((x, y, z))
                else:
                    print(f"Warning: Line '{line.strip()}' is malformed and will be skipped.")

    def adjust_scale(self, target_distance=3.8):
        """
        Calculate a scaling factor based on the average distance between adjacent Cα atoms,
        and adjust to match the target_distance (default 3.8 Å).
        """

        def calculate_average_distance(coords):
            distances = []
            for i in range(len(coords) - 1):
                x1, y1, z1 = coords[i]
                x2, y2, z2 = coords[i + 1]
                distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
                distances.append(distance)
            return np.mean(distances)

        if len(self.coordinates) < 2:
            print("Not enough residues to compute distances.")
            self.scaled_coordinates = self.coordinates[:]
            return

        current_avg_distance = calculate_average_distance(self.coordinates)
        scale_factor = target_distance / current_avg_distance
        self.scaled_coordinates = [(x * scale_factor, y * scale_factor, z * scale_factor)
                                   for x, y, z in self.coordinates]
        print(f"Scale factor applied: {scale_factor:.4f}")

    def write_ca_pdb(self):
        """
        Write the scaled Cα coordinates to a PDB file containing only CA atom information,
        using a standard fixed-width format to ensure MODELLER can parse it.
        """
        res_name_map = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
            'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
            'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
            'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
        }
        with open(self.ca_pdb_filename, 'w') as f:
            # Optional: you可以添加HEADER或MODEL记录以增强兼容性
            # f.write("HEADER    CA-ONLY MODEL\n")
            for i, (res_name, (x, y, z)) in enumerate(zip(self.sequence, self.scaled_coordinates), start=1):
                res_name_3 = res_name_map.get(res_name.upper(), 'UNK')
                # 使用固定宽度格式：参考PDB标准
                line = ("{record:<6s}{serial:>5d} {name:^4s}{alt:1s}{resname:>3s} {chain:1s}"
                        "{resseq:>4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:>2s}\n").format(
                    record="ATOM",
                    serial=i,
                    name="CA",
                    alt="",
                    resname=res_name_3,
                    chain=self.chain_id,
                    resseq=i,
                    x=x, y=y, z=z,
                    element="C")
                f.write(line)
            f.write("END\n")
        print(f"PDB file '{self.ca_pdb_filename}' has been written.")

    def prepare_alignment(self, sequence_name='protein_full'):
        sequence_str = ''.join(self.sequence) + '*'
        start_residue = 1
        end_residue = len(self.sequence)
        with open(self.alignment_filename, 'w') as f:
            # Template section (from our CA file)
            f.write(f">P1;{self.template_code}\n")
            f.write(
                f"structureX:{self.template_code}:{start_residue}:{self.chain_id}:{end_residue}:{self.chain_id}::::\n")
            f.write(sequence_str + "\n\n")
            # Target section (full model)
            f.write(f">P1;{sequence_name}\n")
            f.write(f"sequence:{sequence_name}:{start_residue}:{self.chain_id}:{end_residue}:{self.chain_id}::::\n")
            f.write(sequence_str + "\n")
        print(f"Alignment file '{self.alignment_filename}' has been prepared.")

    def generate_full_model(self):
        log.none()
        env = environ()
        env.io.atom_files_directory = [self.base_dir]

        aln = alignment(env)
        mdl = model(env, file=self.ca_pdb_filename,
                    model_segment=(f'FIRST:{self.chain_id}', f'LAST:{self.chain_id}'))
        aln.append_model(mdl, align_codes=self.template_code, atom_files=self.ca_pdb_filename)
        aln.append(file=self.alignment_filename, align_codes='protein_full')
        aln.align2d()

        class MyModel(automodel):
            def special_patches(self, aln):
                self.rename_segments(segment_ids=['A'])

        a = MyModel(env,
                    alnfile=self.alignment_filename,
                    knowns=self.template_code,
                    sequence='protein_full',
                    assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model = 1
        a.ending_model = 1
        a.make()

        best_model_path = a.outputs[0]['name']
        print("MODELLER built a model:", best_model_path)
        shutil.copyfile(best_model_path, self.full_model_filename)
        print(f"Final full model has been saved as '{self.full_model_filename}'.")

        best_model_path = a.outputs[0]['name']
        print("MODELLER built a model:", best_model_path)
        shutil.copyfile(best_model_path, self.full_model_filename)
        print(f"Final full model has been saved as '{self.full_model_filename}'.")

    # Below are new functions for adding hydrogens, balancing charges, and converting to PDBQT
    def _is_tool_available(self, tool_name):
        """Check if a specified tool (e.g., obabel) is available in the system."""
        return subprocess.call(f"type {tool_name}", shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

    def calculate_center_of_mass(self, structure_file):
        """
        Calculate the center of mass of all atoms in the specified PDB file.
        """
        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure('structure', structure_file)
        except Exception as e:
            raise ValueError(f"Failed to parse file {structure_file}: {e}")
        atom_coords = []
        for atom in structure.get_atoms():
            atom_coords.append(atom.get_coord())
        if not atom_coords:
            raise ValueError(f"No atom coordinates found in file {structure_file}.")
        atom_coords = np.array(atom_coords)
        center_of_mass = atom_coords.mean(axis=0)
        print(f"Center of mass: X={center_of_mass[0]:.3f}, Y={center_of_mass[1]:.3f}, Z={center_of_mass[2]:.3f}")
        return center_of_mass

    def translate_to_origin(self, input_file, output_file):
        """
        Translate the molecule so that its center of mass moves to the origin, then save as a new file.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', input_file)
        center_of_mass = self.calculate_center_of_mass(input_file)
        translation = -center_of_mass  # Calculate the translation vector
        for atom in structure.get_atoms():
            atom.set_coord(atom.get_coord() + translation)
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)
        print(f"Molecule has been translated to the origin and saved as {output_file}")

    def _convert_to_pdbqt_with_obabel(self, input_file, output_pdbqt):
        """Use Open Babel to add hydrogens to the PDB, balance charges, and convert to PDBQT format."""
        command = [
            "obabel", input_file, "-O", output_pdbqt,
            "--partialcharge", "gasteiger", "-h", "-xr"
        ]
        subprocess.run(command, check=True)

    def prepare_pdbqt(self, translate=True, output_translated_file=None):
        """
        Prepare a PDBQT file:
          - Check if Open Babel is available;
          - (Optional) Translate the structure to the origin;
          - Use Open Babel to add hydrogens, compute Gasteiger charges, and generate the PDBQT file.
        """
        # if not self._is_tool_available("obabel"):
        #     raise EnvironmentError("Open Babel is not installed or not in PATH.")

        if translate:
            if output_translated_file is None:
                output_translated_file = self.full_model_filename.replace(".pdb", "_translated.pdb")
            print("Translating molecule to the origin...")
            self.translate_to_origin(self.full_model_filename, output_translated_file)
            input_for_conversion = output_translated_file
        else:
            input_for_conversion = self.full_model_filename

        print("Converting to PDBQT using Open Babel...")
        self._convert_to_pdbqt_with_obabel(input_for_conversion, self.pdbqt_filename)
        print(f"Conversion complete. Intermediate PDBQT file: {self.pdbqt_filename}")

        # If a different docking_folder is specified, move the pdbqt there
        if os.path.abspath(self.docking_folder) != os.path.abspath(self.base_dir):
            final_pdbqt_path = os.path.join(self.docking_folder, os.path.basename(self.pdbqt_filename))
            shutil.move(self.pdbqt_filename, final_pdbqt_path)
            self.pdbqt_filename = final_pdbqt_path
            print(f"PDBQT file has been moved to {self.pdbqt_filename}")

    def prepare_ligand_pdbqt(self, ligand_mol2_path, output_folder=None):
        """
        Process the specified MOL2 ligand file:
          - Read the MOL2 file and calculate its geometric center,
          - Translate all atoms so the geometric center is at the origin,
          - Use Open Babel to convert the translated MOL2 file to PDBQT,
          - Save the resulting ligand PDBQT file to the specified directory (defaults to docking_folder).

        Returns:
        --------
        ligand_pdbqt (str) : Path to the generated ligand PDBQT file
        """
        if output_folder is None:
            output_folder = self.docking_folder

        base_name = os.path.basename(ligand_mol2_path)
        name_no_ext, ext = os.path.splitext(base_name)
        # Temporary path for the translated MOL2 file
        translated_mol2 = os.path.join(self.base_dir, f"{name_no_ext}_translated.mol2")
        # Initial path for the converted ligand PDBQT file in base_dir
        ligand_pdbqt = os.path.join(self.base_dir, f"{name_no_ext}.pdbqt")

        # Use Mol2Translator for translation
        translator = Mol2Translator(ligand_mol2_path, translated_mol2)
        translator.prepare_translated_mol2()

        # Check if Open Babel is available
        # if not self._is_tool_available("obabel"):
        #     raise EnvironmentError("Open Babel is not installed or not in PATH.")

        print("Converting ligand MOL2 to PDBQT using Open Babel...")
        command = [
            "obabel", translated_mol2, "-O", ligand_pdbqt,
            "--partialcharge", "gasteiger", "-h"]
        subprocess.run(command, check=True)
        print(f"Ligand conversion complete. Output pdbqt file: {ligand_pdbqt}")

        # Move the result to output_folder if needed
        if os.path.abspath(output_folder) != os.path.abspath(self.base_dir):
            final_path = os.path.join(output_folder, os.path.basename(ligand_pdbqt))
            shutil.move(ligand_pdbqt, final_path)
            ligand_pdbqt = final_path
            print(f"Ligand PDBQT file has been moved to {ligand_pdbqt}")

        return ligand_pdbqt

    def run_pipeline(self):
        """
        Run the entire pipeline from XYZ to PDBQT in one go (for proteins).
        """
        print("==== 1) Reading XYZ ====")
        self.read_xyz()
        print("==== 2) Adjusting scale ====")
        self.adjust_scale()
        print("==== 3) Writing CA PDB ====")
        self.write_ca_pdb()
        print("==== 4) Preparing alignment ====")
        self.prepare_alignment(sequence_name='protein_full')
        print("==== 5) Generating full model ====")
        self.generate_full_model()
        print("==== 6) Preparing PDBQT file ====")
        self.prepare_pdbqt(translate=True)
        print("==== Pipeline DONE ====")
        print(f"Final docking protein file: {self.pdbqt_filename}")


class Mol2Translator:
    """
    A class to translate the geometric center of a molecule in a MOL2 file to the origin.
    """

    def __init__(self, input_mol2, output_mol2):
        """
        Initialize the Mol2Translator class.

        Parameters:
        - input_mol2 (str): Path to the input MOL2 file.
        - output_mol2 (str): Path to save the translated MOL2 file.
        """
        self.input_mol2 = input_mol2
        self.output_mol2 = output_mol2
        self.atoms = []  # List to store atom information

    def parse_mol2(self):
        """
        Parse the MOL2 file and extract atom information.
        """
        if not os.path.exists(self.input_mol2):
            raise FileNotFoundError(f"The file {self.input_mol2} does not exist.")

        with open(self.input_mol2, 'r') as file:
            lines = file.readlines()

        atom_section = False
        for line in lines:
            line = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                continue
            elif line.startswith("@<TRIPOS>") and atom_section:
                # End of ATOM section
                break

            if atom_section:
                parts = line.split()
                if len(parts) < 6:
                    continue  # Skip malformed lines
                atom_id = parts[0]
                atom_name = parts[1]
                try:
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                except ValueError:
                    raise ValueError(f"Invalid coordinates for atom {atom_id} in file {self.input_mol2}.")
                atom_type = parts[5]
                additional = parts[6:] if len(parts) > 6 else []
                self.atoms.append({
                    'atom_id': atom_id,
                    'atom_name': atom_name,
                    'x': x,
                    'y': y,
                    'z': z,
                    'atom_type': atom_type,
                    'additional': additional
                })

        if not self.atoms:
            raise ValueError(f"No atoms found in the MOL2 file {self.input_mol2}.")

    def calculate_geometric_center(self):
        """
        Calculate the geometric center of the molecule.

        Returns:
        - tuple: (center_x, center_y, center_z)
        """
        sum_x = sum(atom['x'] for atom in self.atoms)
        sum_y = sum(atom['y'] for atom in self.atoms)
        sum_z = sum(atom['z'] for atom in self.atoms)
        num_atoms = len(self.atoms)

        center_x = sum_x / num_atoms
        center_y = sum_y / num_atoms
        center_z = sum_z / num_atoms

        print(f"Geometric center: X={center_x:.3f}, Y={center_y:.3f}, Z={center_z:.3f}")
        return (center_x, center_y, center_z)

    def translate_atoms(self, center):
        """
        Translate all atoms so that the geometric center is at the origin.

        Parameters:
        - center (tuple): The geometric center coordinates (center_x, center_y, center_z).
        """
        center_x, center_y, center_z = center
        for atom in self.atoms:
            atom['x'] -= center_x
            atom['y'] -= center_y
            atom['z'] -= center_z

    def write_translated_mol2(self):
        """
        Write the translated atoms to a new MOL2 file, preserving the original structure.
        """
        with open(self.input_mol2, 'r') as file:
            lines = file.readlines()

        with open(self.output_mol2, 'w') as file:
            atom_section = False
            for line in lines:
                stripped_line = line.strip()
                if stripped_line.startswith("@<TRIPOS>ATOM"):
                    atom_section = True
                    file.write(line)
                    continue
                elif stripped_line.startswith("@<TRIPOS>") and atom_section:
                    atom_section = False
                    file.write(line)
                    continue

                if atom_section:
                    parts = line.split()
                    if len(parts) < 6:
                        file.write(line)
                        continue  # Write malformed lines as is
                    atom_id = parts[0]
                    # Find the corresponding atom in self.atoms
                    atom = next((a for a in self.atoms if a['atom_id'] == atom_id), None)
                    if atom is None:
                        file.write(line)
                        continue
                    new_line = f"{atom['atom_id']:>6} {atom['atom_name']:<10} {atom['x']:>8.3f} {atom['y']:>8.3f} {atom['z']:>8.3f} {atom['atom_type']}"
                    if atom['additional']:
                        additional = ' ' + ' '.join(atom['additional'])
                        new_line += additional
                    new_line += '\n'
                    file.write(new_line)
                else:
                    file.write(line)

        print(f"Translated MOL2 file saved as {self.output_mol2}.")

    def prepare_translated_mol2(self):
        """
        Perform the full preparation: parse, calculate geometric center, translate, and write new MOL2.
        """
        self.parse_mol2()
        center = self.calculate_geometric_center()
        self.translate_atoms(center)
        self.write_translated_mol2()