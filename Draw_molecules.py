# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Dec 23 12:30:41 2024
"""

import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw

# Define the SMILES string for the molecule
smiles = 'CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)CC2=CC(=C(C(=C2)C(C)(C)C)O)C(C)(C)C'

# Convert SMILES to a molecule object
mol = Chem.MolFromSmiles(smiles)

# Draw the molecule as an image
img = Draw.MolToImage(mol, size=(800, 800))  # Increase size for higher resolution

# Display the image using matplotlib
plt.figure(figsize=(6, 6), dpi=350)  # Adjust DPI for high resolution
plt.imshow(img)
plt.axis('off')  # Turn off axis for a cleaner display
plt.tight_layout()
plt.show()

