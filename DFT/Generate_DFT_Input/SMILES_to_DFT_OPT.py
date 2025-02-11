# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Nov 27 01:40:48 2024
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import os

def smiles_to_dft_input(smiles,
                        method="B3LYP",
                        basis_set="6-31G",
                        charge=0,
                        multiplicity=1,
                        filename="dft_opt_input.gjf"):
    
    # Ensure the directory for the file exists
    directory = os.path.dirname(filename)
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    else:
        print('Directory already exist. Skipping!')
        return
        
    # Step 1: Generate a 3D structure
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # Generate 3D coordinates
    AllChem.MMFFOptimizeMolecule(mol)  # Pre-optimize using MMFF
    
    # Step 2: Extract atomic symbols and coordinates
    conf = mol.GetConformer()
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    
    # Step 3: Format Gaussian input
    lines = [
        f"# {method}/{basis_set} Opt",
        "",
        f"Molecule generated from SMILES: {smiles}",
        "",
        f"{charge} {multiplicity}",
    ]
    for atom, coord in zip(atoms, coords):
        lines.append(f"{atom} {coord.x:.6f} {coord.y:.6f} {coord.z:.6f}")
    lines.append("\n")
    
    # Step 4: Write to a file
    with open(filename, "w") as f:
        f.write("\n".join(lines))
    
    print(f"Gaussian input file generated: {filename}")



# User Input
name = "A0001"
smiles = "CCOC(=O)CCC1=CC(=C(O)C(=C1)C(C)(C)C)C(C)(C)C"
method="B3LYP"
basis_set="6-31G"
charge=0
multiplicity=1

dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT'
filename=dirr+rf"\{name}\opt.inp"
smiles_to_dft_input(smiles, method=method, basis_set=basis_set,
                    charge=charge, multiplicity=multiplicity,
                    filename=filename)
