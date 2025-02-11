# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Oct  8 11:32:47 2024
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import magnolia.plot_template as mplt

### Directory ###
AO = 'A0001'
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Production\fixed-300K".format(AO)
filename  = "\\bonds.out"
bondfile = directory+filename
atominfo = bfp.parsebondfile(bondfile, cutoff=0.3, mtypes=True)
#%%

def scav_onset_old(atominfo):
    neighbours = atominfo['neighbours']
    atypes  = atominfo['atypes']
    mtypes  = atominfo['mtypes']

    checklist=set()
    counts=[]
    for step, neigh in neighbours.items():
        molecules=bfp.get_molecules(neigh)
        count=0
        for molecule in molecules:
            moltype=mtypes[list(molecule)[0]]
            if len(molecule)==93 and moltype not in checklist:
                count+=1
                # checklist.add(moltype)
        counts.append(count)
    return counts

def scav_onset(atominfo):
    neighbours = atominfo['neighbours']
    atypes  = atominfo['atypes']
    mtypes  = atominfo['mtypes']
    
    checklist=[]
    for step, neigh in neighbours.items():
        for parent, children in neigh.items():
            if atypes[parent]==3 and mtypes[parent]<=20 and parent not in checklist:
                for child in children:
                    if atypes[child]==1:
                        print(step, parent, 'YES', child)
                        checklist.append(parent)
                        
counts=np.array(scav_onset_old(atominfo))
#%%
steps = np.array(list(atominfo['neighbours'].keys()))/4000

fig, ax = plt.subplots(dpi=350)
plt.rcParams['font.size'] = 15

x = steps.copy()*4+300
y = counts.copy()

# Scatter plot
ax.scatter(x, y, s=1)

# Define the polynomial and logistic functions
def polynomial(x, a, b, c):
    return a * x**2 + b * x + c

# Fit the polynomial function to the data
poly_params, _ = curve_fit(polynomial, x, y, maxfev=10000)

# Generate fits
x_fit = np.linspace(min(x), max(x), 500)
y_poly_fit = polynomial(x_fit, *poly_params)

# Plot data and polynomial fit using ax
ax.scatter(x, y, label='Data', color='blue', s=0.1)
ax.plot(x_fit, y_poly_fit, label='Polynomial fit', color='green')

# Add legend, labels, and title
ax.legend()
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Number of scavenged radicals')

# Solve the quadratic equation for y=5
a, b, c = poly_params
coeffs = [a, b, c - 10]

# Calculate the roots and find the smallest one
onset = np.roots(coeffs).min()

# Plot the vertical line up to y=5 at the step value where y=5
# ax.plot([onset, onset], [0, 10], color='red', linestyle='--', label=f'Step at y=5: {onset:.2f}')
# ax.set_ylim(0, 20)
ax.set_title(AO+"\n")
ax.set_xlim(0)

# Display the plot
plt.show()
print(AO, onset)

#%%
# AOs=['A0001', 'A0002', 'A0003', 'A0004', 'A0005']
# onsets=np.array([75.530830721018, 115.37436618112142, 46.22048068357873,
#         72.31647997010752, 29.751146380216284])*4+300

# fig, ax = plt.subplots(dpi=350)
# ax.bar(AOs,onsets, color='green')
# ax.set_xlabel('Antioxidants')
# ax.set_ylabel('Onset of scavenging (K)')