# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 15 00:51:36 2024
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define constants
timestep = 0.25  # fs
fs_to_ps = 0.1  # 1 fs = 0.1 ps

# Define antioxidants and simulations
AOs = [f'A{i+1:04}' for i in range(5)]
msd_data = []  # List to store MSD data for each AO and each simulation

# Load MSD data for each AO and each simulation
for AO in AOs:
    for sim in range(1, 4):
        dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\{AO}\Equilibration\Sim-{sim}"
        filename = r"\msd_AO_all.txt"
        msd_datafile = dirr + filename

        # Load MSD data
        data = np.loadtxt(msd_datafile, skiprows=1)
        msd_total = data[:, 3]
        
        # Time array in ps
        time = np.array(range(len(msd_total))) * timestep * fs_to_ps
        
        # Add data to list for DataFrame creation
        for t, msd in zip(time, msd_total):
            msd_data.append({"AO": AO, "Time (ps)": t, "MSD (Å^2)": msd, "Simulation": sim})

# Convert to DataFrame for Seaborn
msd_df = pd.DataFrame(msd_data)

# Plot with Seaborn
plt.figure(dpi=350)
sns.lineplot(data=msd_df, x="Time (ps)", y="MSD (Å^2)", hue="AO",
             errorbar=('ci', 95), err_style="band")
plt.ylabel("MSD ($Å^2$)")
plt.xlabel("Time (ps)")
plt.title("Mean Square Displacement (MSD) for Each Antioxidant")
plt.show()

