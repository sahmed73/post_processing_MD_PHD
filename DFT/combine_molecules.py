# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Dec  3 16:02:51 2024
"""

import os
import numpy as np
from ase import Atoms
from ase.io import read, write

# Define the file paths
dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\Try_From_Scratch\Single_Molecules"
molecule1_file = dirr+r"\PAO_Radical.xyz"
molecule2_file = dirr+r"\BHT.xyz"
output_file = dirr+r'\Pre_Reaction.xyz'

# Read molecules from XYZ files
molecule1 = read(molecule1_file)
molecule2 = read(molecule2_file)

# Calculate the bounding box of molecule1
positions1 = molecule1.get_positions()
min1 = positions1.min(axis=0)  # Minimum x, y, z
max1 = positions1.max(axis=0)  # Maximum x, y, z
extent1 = max1 - min1          # Size of molecule1

# Calculate a safe offset for molecule2
buffer = np.array([-12.0]*3)
offset = extent1 + buffer  # Add a 3 Å buffer to avoid overlap
molecule2.translate(offset)

# Combine the two molecules
combined = molecule1 + molecule2

# Save the combined structure to a new XYZ file
write(output_file, combined)

print(f"Combined molecule written to {output_file}")
