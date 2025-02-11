# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Apr  7 14:42:35 2024
"""

import pandas as pd
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt
mplt.custom_plot_features()

# Adjust the path based on where you've stored your files in Google Drive
path_to_files = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\OPLSAA\PAO_Antioxidant\MSD'

# Load the data
df_a = pd.read_csv(path_to_files + '\\A_msd_antioxidant_data.txt', sep=r"\s+")
df_b = pd.read_csv(path_to_files + '\\B_msd_antioxidant_data.txt', sep=r"\s+")
df_bht = pd.read_csv(path_to_files + '\\BHT_msd_antioxidant_data.txt', sep=r"\s+")

# Since the time steps are the same, you can just use the time steps from one file
time = df_a['Timestep']
time = time-time.min()
time = time*0.25/1000


plt.plot(time, df_a['MSD'], label='A')
plt.plot(time, df_b['MSD'], label='B')
plt.plot(time, df_bht['MSD'], label='BHT')

plt.xlabel('Time (ps)')
plt.ylabel('MSD')
# plt.title('MSD vs. Time Steps for A, B, and BHT\n')
plt.legend(loc='upper left')
mplt.saveplot(name='MSD_A_B_BHT')
plt.show()