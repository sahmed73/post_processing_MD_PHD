# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  1 14:23:22 2024
"""

import numpy as np
import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt
import pandas as pd
from sklearn.metrics import r2_score

result = {}
for sim in range(1):
    dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0002\Production\Sim-{}".format(1+sim)
    filename=r"\bonds.out"
    bondfile=dirr+filename
    atominfo=bfp.parsebondfile(bondfile,mtypes=True)
    
    neighbours=atominfo["neighbours"]
    atypes=atominfo["atypes"]
    mtypes=atominfo["mtypes"]
    
    H, C, O = 1, 2, 3 # atom types
    timestep = 0.25
    N_PAOr = 20 # number of radicals
    N_AO = 15 # number of antioxidants
    
    SR = {} # count of scavenged radicals
    for step, neigh in neighbours.items():
        molecules=bfp.get_molecules(neigh)
        G = nx.Graph(neigh)
        time = step*timestep/1000 # in ps
        SR[time] = 0
        for molecule in molecules:
            g = G.subgraph(molecule)
            # L = [atypes[x] for x in molecule].count(C) # only counting carbons
            L = 0
            for atom in molecule:
                if atypes[atom]==C and mtypes[atom]<=N_PAOr:
                    L+=1
            for parent in g.nodes():
                if atypes[parent]==O and mtypes[parent]>N_PAOr:
                    SR[time]+=L
    result[f'Sim-{sim+1}']=SR.copy()
#%%
plt.rcParams['font.size']=13
df  = pd.DataFrame(result)
dff = df.dropna(how='any')/N_PAOr
fig, ax =plt.subplots(dpi=350)   
mplt.mean_range_plot(dff, ax=ax)
# ax.scatter(SR.keys(), SR.values())
ax.set_xlabel("Time (ps)")
ax.set_ylabel("Number of SR")

# Define ramp rate and initial temperature
ramp_rate = 4  # Temperature change per unit time
initial_temp = 300  # Starting temperature

fig, ax = plt.subplots(dpi=350)

# Primary plot on the original axis
mplt.mean_range_plot(dff, ax=ax, label='Actual Data',color='tab:orange')
ax.set_xlabel("Time (ps)")
ax.set_ylabel("Number of SR")

# Create a secondary x-axis for temperature
secax = ax.secondary_xaxis('top', functions=(lambda x: x * ramp_rate + initial_temp, 
                                             lambda x: (x - initial_temp) / ramp_rate))
secax.set_xlabel("Temperature (K)")
ax.grid(alpha=0.1)

# ## fitting
# # Extract the x and y values for fitting
# time = dff.index # Replace 'Time (ps)' with the actual column name if different
# number_of_sr = dff.mean(axis=1) # Replace 'Number of SR' with the actual column name if different

# # Perform a quadratic fit
# coefficients = np.polyfit(time, number_of_sr, 2)  # 2 represents the degree of the polynomial (quadratic)
# quadratic_fit = np.poly1d(coefficients)

# # Generate x values for plotting the fit
# x_fit = np.linspace(time.min(), time.max(), 100)
# y_fit = quadratic_fit(x_fit)

# # Plot the quadratic fit
# ax.plot(x_fit, y_fit, color='red', linestyle='--', label='Quadratic Fit')
# ax.legend(loc='upper right')
# ax.set_ylim(0-1,20+1)
# ax.set_xlim(0-10,175+10)

# # Calculate y value for x = 50 using the quadratic fit
# x_value = 50
# y_value = np.round(quadratic_fit(x_value))
# print(f"The y value at x = {x_value} is {y_value:.2f}")
# plt.plot([x_value,x_value], [0-1,y_value], '--', color='k', alpha=0.5)
# plt.plot([-10,x_value], [y_value,y_value], '--', color='k', alpha=0.5)
# ax.plot(x_value, y_value, 'o', color='k', alpha=0.5)
# ax.text(x_value/2, y_value, f"{y_value:.1f}", ha='right',
#         va='bottom', color='k', alpha=1.0)

# # R-squared
# y_fit_values = quadratic_fit(time)
# r_squared = r2_score(number_of_sr, y_fit_values)
# # Display R^2 value in the plot area
# ax.text(0.05, 0.95, f"$R^2 = {r_squared:.3f}$", transform=ax.transAxes, 
#         verticalalignment='top', horizontalalignment='left', 
#         bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
