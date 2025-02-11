# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Aug 26 18:26:48 2024
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles

# Define the SMILES string
smiles = 'CCCCCCCCCCC(O)(CCCCCCCC)CC(C)CCCCCCCC'

# Convert SMILES to a molecule object
mol = Chem.MolFromSmiles(smiles)

# Add hydrogens to the molecule
mol_with_h = Chem.AddHs(mol)

# Generate 3D coordinates
AllChem.EmbedMolecule(mol_with_h)
AllChem.UFFOptimizeMolecule(mol_with_h)

# Write the molecule to a PDB file
pdb_filename = 'molecule.pdb'
rdmolfiles.MolToPDBFile(mol_with_h, pdb_filename)

print(f"PDB file saved as {pdb_filename}")
