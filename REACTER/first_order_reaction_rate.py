# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Feb 24 11:10:48 2025
"""

import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit

plt.style.use('default')

parent_dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001'

AOs = ['A0001', 'A0002', 'A0003', 'A0004', 'A0005']
MAP = dict(zip(AOs, ['L1','L2','S1','S2','S3']))
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

# First-order reaction function
def first_order_growth(t, k):
    N0 = 50 # initial number of antiocidant
    return N0 * (1 - np.exp(-k * t))

data = {}
for AO in AOs:
    data[AO] = {}
    for i in range(3):  # 3 simulations per antioxidant
        dirr = parent_dirr + rf'\{AO}\Reaction\TSD=1.95\Sim-{i+1}'
        filename = r'\log.lammps'
        logfile = dirr + filename
        thermo = lfp.thermo_panda(logfile, serial=1)

        # Extract time and number of scavenged radicals
        x = thermo['Time'] / 1000  # Convert to ns
        y = thermo['v_rxn1']
        
        if AO in ['A0002', 'A0005']:
            y += thermo['v_rxn2']
        #     # y = y/2
            
        # Fit the function
        try:
            popt, _ = curve_fit(first_order_growth, x,y)  # Increased function evaluations
        
        except RuntimeError:
            print(f"Curve fitting failed for {AO}, Sim-{i+1}. Trying different initial values...")
            try:
                popt, _ = curve_fit(first_order_growth, x, y, 
                                     p0=[0.05, max(y) * 0.8], bounds=(0, np.inf), 
                                     maxfev=20000)
            except RuntimeError:
                print(f"Final attempt failed for {AO}, Sim-{i+1}. Skipping this dataset.")
                popt = [np.nan, np.nan]  # Assign NaN if fitting fails

        # Extract fitted parameters
        k_fit = popt[0]
        print(f"Fitted rate constant (k): {k_fit:.4f}")

        data[AO][f'Sim-{i+1}'] = k_fit  # Store the reaction rate constant
#%%
plt.style.use('default')
plt.rcParams['font.size']=22

# Convert to DataFrame
df = pd.DataFrame(data).T

fig, ax = plt.subplots(dpi=350)
x = [MAP[AO] for AO in df.index]  # X-axis: Antioxidant names
y = df.mean(axis=1)       # Y-axis: Mean k values
err = df.std(axis=1) / np.sqrt(df.shape[1])  # n = number of simulations (3)


ax.bar(x,y,yerr=err,capsize=5,color=colors)
ax.set_xlabel('Antioxidants')
ax.set_ylabel('$K_{Scavenging}$ (ns$^{-1}$)')
