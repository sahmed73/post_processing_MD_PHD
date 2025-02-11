# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 28 03:28:04 2024
"""

import numpy as np
import matplotlib.pyplot as plt

dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Equilibration\Sim-1"
filename=r"\msd_AO_all.txt"
msd_datafile=dirr+filename

# Load the data from the text file (assuming it's named 'msd_data.txt')
data = np.loadtxt(msd_datafile, skiprows=1)

# Split the data into separate variables
msd_x = data[:, 0]
msd_y = data[:, 1]
msd_z = data[:, 2]
msd_total = data[:, 3]

# Plotting the MSD components
plt.figure(dpi=120)
plt.plot(msd_x, label='MSD X')
plt.plot(msd_y, label='MSD Y')
plt.plot(msd_z, label='MSD Z')
plt.plot(msd_total, label='MSD Total', linestyle='--', linewidth=2)

# Adding labels and legend
plt.xlabel('Time Step')
plt.ylabel('MSD')
plt.title('Mean Squared Displacement (MSD) vs. Time')
plt.legend()

# Show the plot
plt.show()

