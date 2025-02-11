# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  1 12:54:53 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Define constants
timestep = 0.25  # fs
fs_to_ps = 0.1  # 1 fs = 0.1 ps

# Define antioxidants and simulations
AOs = ['A0001', 'A0002', 'A0003']
num_simulations = 8
diffusivities = {AO: [] for AO in AOs}  # Dictionary to store diffusivities for each AO

# Calculate diffusivity for each AO and each simulation
for AO in AOs:
    for sim in ['Sim-1','NVT-2','NVT-3']:
        dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\{AO}\Equilibration\{sim}"
        filename = r"\msd_AO_all.txt"
        msd_datafile = dirr + filename

        # Load MSD data
        data = np.loadtxt(msd_datafile, skiprows=1)
        msd_total = data[:, 3]
        
        # Time array in ps
        time = np.array(range(len(msd_total))) * timestep * fs_to_ps
        
        # Calculate diffusivity using linear regression
        slope, intercept, r_value, p_value, std_err = linregress(time, msd_total)
        diffusivity = slope / 6  # Einstein relation for 3D diffusion
        diffusivities[AO].append(diffusivity)  # Store diffusivity for each simulation

        print(f"Diffusivity for {AO}, Sim-{sim}: {diffusivity:.4e} Å²/ps")

# Calculate average diffusivity and standard deviation for each AO
avg_diffusivities = [np.mean(diffusivities[AO]) for AO in AOs]
std_diffusivities = [np.std(diffusivities[AO]) for AO in AOs]

# Plot diffusivity bar plot
fig, ax = plt.subplots(dpi=350)
ax.bar(AOs, avg_diffusivities, yerr=std_diffusivities, color=plt.cm.tab10.colors[:len(AOs)], capsize=5, width=0.7)
ax.set_ylabel("Diffusivity ($Å^2/ps$)")
ax.set_title("Average Diffusivity for Each AO with 8 Simulations")

plt.show()
