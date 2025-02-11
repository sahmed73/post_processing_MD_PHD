# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  8 04:26:37 2024
"""


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Define constants
timestep = 0.25  # fs
fs_to_ps = 0.1  # 1 fs = 0.1 ps

# Define antioxidants and simulations
AOs = [f'A{i+1:04}' for i in range(5)]
diffusivities_data = []  # List to store diffusivities and associated metadata for each AO

# Calculate diffusivity for each AO and each simulation
for AO in AOs:
    for sim in range(1,4):
        dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\{AO}\Equilibration\Sim-{sim}"
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
        diffusivities_data.append({"AO": AO, "Diffusivity": diffusivity, "Simulation": sim})
        print(f"Diffusivity for {AO}, {sim}: {diffusivity:.4e} Å²/ps")

# Convert to DataFrame for Seaborn
diffusivities_df = pd.DataFrame(diffusivities_data)

# Plot with Seaborn
plt.figure(dpi=350)
sns.barplot(data=diffusivities_df, x="AO", y="Diffusivity", ci=95, capsize=0.1, palette="tab10")
plt.ylabel("Diffusivity ($Å^2/ps$)")
plt.xlabel('')
# plt.title("Average Diffusivity for Each AO with Confidence Intervals")
plt.show()
