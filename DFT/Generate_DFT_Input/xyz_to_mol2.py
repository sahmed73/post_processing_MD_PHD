# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Dec  3 14:31:27 2024
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Define the file paths
dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\A0001\pre'
xyzfile = os.path.join(dirr, 'opt_input.out.xyz')
mol2file = xyzfile.replace('.xyz', '.mol2')

# Function to read XYZ and convert to RDKit molecule
def xyz_to_rdkit_mol(xyzfile):
    with open(xyzfile, 'r') as f:
        lines = f.readlines()

    # Parse the XYZ file
    atom_count = int(lines[0].strip())
    atoms_data = lines[2:2 + atom_count]

    mol = Chem.RWMol()
    coords = []

    for line in atoms_data:
        element, x, y, z = line.split()
        atom = Chem.Atom(element)
        idx = mol.AddAtom(atom)
        coords.append((float(x), float(y), float(z)))

    conf = Chem.Conformer(len(coords))
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf)

    return mol

# Convert XYZ to RDKit molecule
mol = xyz_to_rdkit_mol(xyzfile)

# Add bond information if needed (optional)
AllChem.EmbedMolecule(mol)

# Write the RDKit molecule to a MOL2 file
Chem.MolToMol2File(mol, mol2file)
print(f"Converted {xyzfile} to {mol2file}")


