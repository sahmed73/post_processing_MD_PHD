# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  8 13:58:55 2024
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Load the Rg data from CSV
plt.rcParams['font.size'] = 16
fig, ax = plt.subplots(dpi=350)

AOs = [f'A{i+1:04}' for i in range(5)]
for AO in AOs:
    dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\python_codes"
    rg_data = pd.read_csv(dirr + f"\\rg_values_{AO}.csv")
    
    # Plot histogram without KDE
    # sns.histplot(rg_data["Rg"], kde=False, stat="density", bins=30, alpha=0.6, label=AO, ax=ax)
    
    # Plot KDE separately with consistent line width
    sns.kdeplot(rg_data["Rg"], ax=ax, linewidth=1.0, shade=True, label=AO)  # Adjust `linewidth` as needed
    
plt.xlabel("Radius of Gyration (Å)")
plt.ylabel("Density")
plt.title("Probability Density of Radius of Gyration\n")
plt.legend()
plt.show()

