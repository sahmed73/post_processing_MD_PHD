# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Mar 14 10:10:57 2025
"""

import MDAnalysis as mda
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt

dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0001\Eq\Interaction_Energy+MSD"
topofile = dirr+r'\\min.data'
dumpfile = dirr+"\\eq.nvt.unwrapped.lammpstrj"

# Load MD trajectory
u = mda.Universe(topofile, dumpfile, format='LAMMPSDUMP')

# Define antioxidant molecule volume (adjust for your AO type)
mol = Chem.MolFromSmiles("CCOC(=O)CCC1=CC(=C(O)C(=C1)C(C)(C)C)C(C)(C)C")
mol_PAO = Chem.MolFromSmiles("CCCCCCCCCCC([O])(CCCCCCCC)CC(C)CCCCCCCC")
V_vdw = (Descriptors.MolMR(mol)+Descriptors.MolMR(mol_PAO))*1.3  # Convert to Å³ using Bondi's method

# Iterate over all frames
ffv_values = []
for ts in u.trajectory:
    # Extract box dimensions (Total volume V)
    Lx, Ly, Lz = ts.dimensions[:3]
    V_total = Lx * Ly * Lz

    # Count the number of antioxidant molecules
    num_AO = 150 #len(u.select_atoms("resid 101:150"))  # Adjust selection
    # Calculate occupied volume
    V_occupied = num_AO * V_vdw

    # Compute FFV
    FFV = 1 - (V_occupied / V_total)
    ffv_values.append(FFV)

    print(f"Frame {ts.frame}: FFV = {FFV:.4f}")

# Compute average FFV over all frames
avg_FFV = np.mean(ffv_values)
print(f"Average FFV over trajectory: {avg_FFV:.4f}")
#%%
fig, ax = plt.subplots(dpi=350)
ax.plot(ffv_values)
