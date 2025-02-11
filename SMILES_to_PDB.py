# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed May  1 16:44:54 2024
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import sys

def smiles2pdb(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol_with_h = Chem.AddHs(mol)

    # Generate 3D coordinates for the molecule with hydrogens
    AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_with_h)

    # Convert the molecule to a PDB format
    pdb_block = Chem.MolToPDBBlock(mol_with_h)
    
    return pdb_block

### for command line input
name   = sys.argv[1]
smiles = sys.argv[2]

pdb = smiles2pdb(smiles)
# save pdb
with open(f"{name}.pdb", "w") as file:
    file.write(pdb)