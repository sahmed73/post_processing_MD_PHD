# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov  6 22:12:04 2023
"""

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# Create a molecule object from a SMILES string (you can replace this with your molecule)
mol = Chem.MolFromSmiles('CCO')  # Ethanol as an example

# Compute 2D coordinates for the molecule
AllChem.Compute2DCoords(mol)

# Create a list to store text labels for atoms
atom_labels = [atom.GetSymbol() for atom in mol.GetAtoms()]

# Draw the molecule with atom labels
img = Draw.MolToImage(mol, size=(300, 300), wedgeBonds=True, kekulize=True, wedgeFontSize=0.8, atomLabels=atom_labels)

# Display the image
img.show()
