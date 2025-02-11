# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 28 01:43:55 2024
"""

import magnolia.bondfile_parser as bfp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

### Directory ###
AO='A0002'
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\{}\Production\Sim-1".format(AO)
filename  = "\\bonds.out"
bondfile = directory+filename
atominfo = bfp.parsebondfile(bondfile, mtypes=True)
#%%
def count_scavenged(atominfo):
    neighbours = atominfo['neighbours']
    atypes  = atominfo['atypes']
    mtypes = atominfo['mtypes']
    
    # types
    hydrogen_type=1
    carbon_type=2
    oxygen_type=3
    
    scavenged = {}
    
    scavenged_with_AO_hydroxyl_h={}
    scavenged_with_other_h={}
    
    
    # is the H from AO hydroxyl
    is_from_AO_hydroxyl={} # bool
    for step, neigh in neighbours.items():
        for parent, children in neigh.items():
            is_from_AO_hydroxyl[parent]=False
            
            if atypes[parent]==hydrogen_type and mtypes[parent]>init_num_of_radicals:
                child = children[0]
                if atypes[child]==oxygen_type:
                    is_from_AO_hydroxyl[parent]=True
                if len(children)>1:
                    print("Warning: some H have more than one child")
        break # only checking in the first step
                
    for step, neigh in neighbours.items():
        scavenged[step]=0
        scavenged_with_AO_hydroxyl_h[step]=0
        scavenged_with_other_h[step]=0
        
        for parent, children in neigh.items():
            molecules = bfp.get_molecules(neigh)
            # count if it lookes like -C-O-H
            if atypes[parent]==oxygen_type and mtypes[parent]<=init_num_of_radicals:
                children_types=sorted([atypes[x] for x in children])
                if children_types==[hydrogen_type,carbon_type]: # H and C
                    for molecule in molecules:
                        if parent in molecule:
                            if len(molecule)>=parmissible:
                                scavenged[step]+=1
                            break
                    
                    for child in children:
                        if atypes[child]==hydrogen_type:
                            if is_from_AO_hydroxyl[child] is True:
                                scavenged_with_AO_hydroxyl_h[step]+=1
                            else:
                                scavenged_with_other_h[step]+=1
    
    return pd.Series(scavenged), pd.Series(scavenged_with_AO_hydroxyl_h), pd.Series(scavenged_with_other_h)


timestep  = 0.25
ramp_rate = 4
initial_temp = 300
init_num_of_radicals = 20
init_num_of_AOs = 15
parmissible  = 0
 
scavenged, hydroxyl, others =count_scavenged(atominfo)
time=scavenged.index*timestep/1000
temp=initial_temp+ramp_rate*time

scavenged.index=temp
hydroxyl.index=temp
others.index=temp

plt.rcParams['font.size']=16
fig, ax = plt.subplots(dpi=350)

# Quadratic fit using the index as x-values
x,y = temp, scavenged.values
coefficients = np.polyfit(x, y, deg=2)  # Fit a quadratic polynomial
quadratic_fit = np.poly1d(coefficients)  # Create a polynomial function
y_fit = quadratic_fit(x)

# R-Square
ss_res = np.sum((y - y_fit) ** 2)  # Residual sum of squares
ss_tot = np.sum((y - np.mean(y)) ** 2)  # Total sum of squares
r_squared = 1 - (ss_res / ss_tot)
# Add R^2 text to the plot
ax.text(0.05, 0.95, f'R² = {r_squared:.2f}', transform=plt.gca().transAxes, 
         fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))

ax.scatter(x[::5], y[::5], label='Original Data', color='tab:orange', alpha=0.9, s=30,
           edgecolor='k', linewidths=0.3)
ax.plot(x, y_fit, label='Quadratic Fit', color='red')
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Number of scavenged radicals")
ax.legend(loc='lower right', fontsize=12)
ax.set_title(AO)

# onset
target_radicals = 10
x_fine = np.linspace(min(x), max(x), 1000)
y_fine = quadratic_fit(x_fine)
onset  = np.round(x_fine[np.where(y_fine >= target_radicals)[0][0]])
ax.plot([onset,onset],[-2,target_radicals],'--', color='k', alpha=0.5)
ax.plot([215 ,onset],[target_radicals,target_radicals],'--', color='k', alpha=0.5)
ax.set_ylim(-2,20+2)
ax.set_xlim(215)

ax.text(onset, target_radicals / 2, f'{onset:.1f}', fontsize=12, color='k', 
        verticalalignment='top', horizontalalignment='left', rotation=0)
#%%