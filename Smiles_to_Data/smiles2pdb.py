# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon May 27 03:49:24 2024
"""

import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def main():
    # inputs from the bash script
    inputfile  = sys.argv[1] # smiles.in dirr
    
    smiles = {}
    if os.path.isfile(inputfile):
        with open(inputfile, 'r') as ipf:
            # Skip the first two lines
            next(ipf)
            next(ipf)
            
            for i, line in enumerate(ipf):
                if i % 2 == 0:
                    species = line.strip()
                else:
                    smiles[species] = line.strip()
    else:
        for i, item in enumerate(sys.argv):
            if i==0: continue
        
            if i % 2 == 1:
                species = sys.argv[i]
            else:
                smiles[species] = sys.argv[i]
        
    ## smiles to single pdb
    for species, smile in smiles.items():
        mol = Chem.MolFromSmiles(smile)
        mol_with_h = Chem.AddHs(mol)
    
        # Generate 3D coordinates for the molecule with hydrogens
        AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol_with_h)
        
        # Convert the molecule to a PDB format
        pdb_block = Chem.MolToPDBBlock(mol_with_h)
        
        with open(f"//{species}.pdb", "w") as file:
            file.write(pdb_block)
            
if __name__=="__main__":
     main()
