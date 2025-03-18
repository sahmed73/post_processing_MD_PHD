# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 28 03:33:50 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Constants
timestep = 1.0  # fs
fs_to_ps = 0.1   # 1 fs = 0.1 ps

AOs = ['A0001', 'A0002', 'A0003', 'A0004', 'A0005']
folders = ["Sim-4_Interaction_Energy+MSD",
           "Sim-1_Interaction_Energy+MSD",
           "Sim-10_Interaction_Energy+MSD",
           "Sim-1_Interaction_Energy+MSD",
           "Sim-1_Interaction_Energy+MSD"]

diffusivities = []

fig, ax = plt.subplots(dpi=350)

for AO, folder in zip(AOs, folders):
    dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\{AO}\Eq\{folder}"
    filename = r"\msd_AO_all.txt"
    msd_datafile = dirr + filename

    # Load data
    data = np.loadtxt(msd_datafile, skiprows=1)
    msd_total = data[:, 3]

    # Convert time to ps
    time = np.array(range(len(msd_total))) * timestep * fs_to_ps

    # Plot full MSD curve
    ax.plot(time, msd_total, label=AO, linestyle='--', linewidth=2)

    # Select the **linear region** of MSD
    min_time = 100  # Lower bound of the linear region (adjust as needed)
    max_time = 500  # Upper bound of the linear region (adjust as needed)
    
    linear_indices = np.where((time >= min_time) & (time <= max_time))[0]
    
    time_linear = time[linear_indices]
    msd_linear = msd_total[linear_indices]

    # Linear regression on selected region
    slope, intercept, r_value, p_value, std_err = linregress(time_linear, msd_linear)
    diffusivity = slope / 6  # Einstein relation for 3D diffusion
    diffusivities.append(diffusivity)

    print(f"Diffusivity for {AO}: {diffusivity:.4e} Å²/ps (Using {min_time}ps to {max_time}ps)")

# Add labels and legend
plt.xlabel('Time (ps)')
plt.ylabel('MSD ($Å^2$)')
plt.legend()
plt.show()

# Plot bar chart for diffusivities
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
fig, ax = plt.subplots(dpi=350)
ax.bar(AOs, diffusivities, width=0.7, color=colors)
ax.set_ylabel("Diffusivity ($Å^2/ps$)")
