# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Nov 27 03:07:49 2024

-------------------------------------------------------------------------------
SMILES Input tamplate
-------------------------------------------------------------------------------
Date 7/11/2024 # comment line 1
smiles inputs  # comment line 2
A0001
CCOC(=O)CCC1=CC(=C(O)C(=C1)C(C)(C)C)C(C)(C)C
A0002
CC(C)(C)C1=CC(\C=C/CO)=CC(=C1O)C1=CC(\C=C/CO)=CC(=C1O)C(C)(C)C
A0003
Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1
A0004
CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)C=CC(=O)
A0005
CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)CC2=CC(=C(C(=C2)C(C)(C)C)O)C(C)(C)C
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
SMILES Input 2.0 tamplate
-------------------------------------------------------------------------------
SMILES # comment line 1
Date 11/17/2024 # comment line 2
# Set-1: PAO Radical
PAO-OH 20
CCCCCCCCCCC(O)(CCCCCCCC)CC(C)CCCCCCCC
# Set-2: Antioxidants
A0001 15
CCOC(=O)CCC1=CC(=C(O)C(=C1)C(C)(C)C)C(C)(C)C
A0002 15
CC(C)(C)C1=CC(\C=C/CO)=CC(=C1O)C1=CC(\C=C/CO)=CC(=C1O)C(C)(C)C
A0003 15
Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1
A0004 15
CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)C=CC(=O)
A0005 15
CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)CC2=CC(=C(C(=C2)C(C)(C)C)O)C(C)(C)C
-------------------------------------------------------------------------------
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


def guess_charge_and_multiplicity(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string")
    
    # Calculate formal charge
    charge = Chem.GetFormalCharge(mol)
    
    # Determine multiplicity
    unpaired_electrons = sum([atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()])
    multiplicity = unpaired_electrons + 1  # Multiplicity = 2S + 1
    
    return charge, multiplicity

def read_smiles_input(smiles_infile):
    with open(smiles_infile) as infile:
        data={}
        next(infile) # skipping first line
        next(infile) # skipping second line
        
        for unstrip_line in infile:
            line=unstrip_line.strip()
            if line.startswith('#'): continue
        
            name=line
            smiles=next(infile)
            data[name]=smiles
    
    return data


smile_path=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\Job_Submission_Scripts_and_Inputs\Automation_&_input_eq-3NPT-NVT_Pro-NVT_Slow_Ramp\Automation_Scripts\smiles.in"
data=read_smiles_input(smile_path)

dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT'
for name, smiles in data.items():
    method="B3LYP"
    basis_set="6-31G"
    charge, multiplicity= guess_charge_and_multiplicity(smiles)
    print(charge, multiplicity)
    filename=dirr+rf"\{name}\opt.inp"
    smiles_to_dft_input(smiles, method=method, basis_set=basis_set,
                        charge=charge, multiplicity=multiplicity,
                        filename=filename)