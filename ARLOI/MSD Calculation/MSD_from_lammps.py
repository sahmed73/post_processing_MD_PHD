# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Oct  4 03:23:28 2024
"""

import pandas as pd
import matplotlib.pyplot as plt

# Assume the file "msd_AO_all.txt" is available in the required format
# Simulate reading a file in the format: "msd_x msd_y msd_z msd_total"
file_path = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\Test_Run\test-1\msd_PAO_all.txt"

# Read the file into a dataframe, assuming it is space-separated
df = pd.read_csv(file_path, sep="\s+")

# Add a "time" column based on the index (which represents timesteps)
df['time'] = df.index

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(df['time'], df['msd_x'], label='MSD X')
plt.plot(df['time'], df['msd_y'], label='MSD Y')
plt.plot(df['time'], df['msd_z'], label='MSD Z')
plt.plot(df['time'], df['msd_total'], label='MSD Total', linewidth=2, linestyle='--')

# Add labels and title
plt.xlabel('Timesteps')
plt.ylabel('MSD')
plt.title('MSD vs Time (Timesteps)')
plt.legend()

# Display the plot
plt.grid(True)
plt.show()
