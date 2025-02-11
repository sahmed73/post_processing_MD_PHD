# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Feb 10 15:56:03 2025
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import sys

def smiles_to_xyz(smiles, name, filename):
    # Convert SMILES to RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogen atoms

    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

    # Get atomic symbols and 3D coordinates
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()

    xyz_lines = [f"{num_atoms}\n{smiles}"]  # XYZ header
    for i in range(num_atoms):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        xyz_lines.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

    # Save to XYZ file
    with open(filename, "w") as file:
        file.write("\n".join(xyz_lines) + "\n")

### Command Line Input Handling
if __name__ == "__main__":
    print("Usage: python script.py molecule_name SMILES")

    name = 'BHT'
    smiles = 'Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1'
    
    dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\Try_From_Scratch\Single_Molecules"
    filename = f"\\{name}.xyz"
    smiles_to_xyz(smiles, name, dirr+filename)
    print(f"XYZ file '{name}.xyz' generated successfully!")