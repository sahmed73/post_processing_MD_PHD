# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  8 04:53:53 2024
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# Define the file paths
trej = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0001\Production\300-500K_TRR=1Kpps\Sim-1\prod.nvt.unwrapped.lammpstrj"
topo = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0001\Equilibration\Sim-1\H_removed_eq_1.data"

# Load the trajectory
u = mda.Universe(topo, trej, format="LAMMPSDUMP")

rg_avg_all = []
for res in range(21,36):    
    selection = u.select_atoms(f"byres (resid {res})")
    
    # Initialize list to store Rg values
    rg_values = []
    
    # Calculate Rg over the trajectory
    for ts in u.trajectory:
        rg = selection.radius_of_gyration()
        rg_values.append(rg)
    
    # Compute the average Rg over the trajectory
    average_rg = np.mean(rg_values)
    std_rg = np.std(rg_values)
    rg_avg_all.append(average_rg)
    print(f"Avg Rg for molecules {res}: {average_rg:.3f} Å")

# Plot Rg over time
time = np.array([ts.time for ts in u.trajectory])/4000
fig, ax = plt.subplots(dpi=350)
# ax.plot(time, rg_values, label='Rg over time')
ax.plot(range(21,36), rg_avg_all)
plt.xlabel('Time (ps)')
plt.ylabel('Radius of Gyration (Å)')
# plt.title('Radius of Gyration (Rg) over Time for Selected Molecules (IDs 21-35)')
plt.legend()
plt.show()

