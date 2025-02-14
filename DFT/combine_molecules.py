# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Dec  3 16:02:51 2024
"""

import numpy as np
import os
from ase.io import read, write

# Define the file paths
dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0005\DataFile\Single_Molecules"

number  = 100 # hard coding . look for max 100, more than enough

for number in range(1,100): # number
    for prepost in ['pre', 'post']:
        
        molecule1_file = dirr+rf"\PAOr_optimized_{prepost}{number}.xyz"
        molecule2_file = dirr+rf"\A0005_optimized_{prepost}{number}.xyz"
        output_file = dirr+rf'\{prepost}_reaction{number}.xyz'
        
        if not os.path.exists(molecule1_file) or not os.path.exists(molecule1_file):
            continue
        
        
        # Read molecules from XYZ files
        molecule1 = read(molecule1_file)
        molecule2 = read(molecule2_file)
        
        # Calculate the bounding box of molecule1
        positions1 = molecule1.get_positions()
        min1 = positions1.min(axis=0)  # Minimum x, y, z
        max1 = positions1.max(axis=0)  # Maximum x, y, z
        extent1 = max1 - min1          # Size of molecule1
        
        # Calculate a safe offset for molecule2
        buffer = np.array([-20.0]*3)
        offset = extent1 + buffer  # Add a 3 Å buffer to avoid overlap
        molecule2.translate(offset)
        
        # Combine the two molecules
        combined = molecule1 + molecule2
        
        # Save the combined structure to a new XYZ file
        write(output_file, combined)
        
        print(f"Combined {prepost}-{number}'s")
