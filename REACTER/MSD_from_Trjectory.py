# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Feb 23 21:19:55 2025
"""


import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AOs          = ['A0001', 'A0002', 'A0003', 'A0004', 'A0005']
n_atoms      = [52, 64, 40, 43, 75]
init_atomids = [9211, 9225, 9209, 9210, 9210]
msd_data = pd.DataFrame()
diffusivity_data = {}

for i, AO in enumerate(AOs):    
    topo = rf'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\{AO}\DataFile\Bulk\100_PAOr_50_{AO}_density=0.2.data'
    trej = rf'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\{AO}\Eq\Interaction_Energy+MSD\eq.nvt.unwrapped.lammpstrj'
    u = mda.Universe(topo, trej, format='LAMMPSDUMP')        
    
    # Select all atoms or a specific group of atoms
    ids = [init_atomids[i]+x*n_atoms[i] for x in range(50)]
    selection_string = " or ".join(f"id {i}" for i in ids)
    AO_atoms = u.select_atoms(selection_string)
    # single_atom = u.select_atoms('id ')
     
    # Compute MSD using the Einstein relation
    msd = EinsteinMSD(AO_atoms, msd_type='xyz', fft=True, unwarp=True).run()
    
    msd_timeseries = msd.results.timeseries
    timestep = 1
    time = np.array([ts.frame * timestep for ts in u.trajectory])
    
    msd_data['Time'] = time
    msd_data[AO] = msd_timeseries
    
#%%
plt.style.use('default')
plt.rcParams['font.size']=17
fig, ax = plt.subplots(dpi=350)
for col in msd_data:
    if col=='Time': continue
    ax.plot(time/1000, msd_data[col], label=col, linestyle='--', linewidth=2)

ax.set_xlabel('Time (ns)')
ax.set_ylabel('MSD ($Å^2$)')
ax.legend()
#%%
dfs = []
for col in msd_data:
    if col=='Time': continue
    m, c = np.polyfit(msd_data['Time']/1000, msd_data[col], 1)
    D = m/6
    dfs.append(D)

# Plot bar chart for diffusivities
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
fig, ax = plt.subplots(dpi=350)
ax.bar(AOs, dfs, width=0.7, color=colors, edgecolor='k', linewidth=1)
ax.set_ylabel("Diffusivity ($Å^2/ns$)")
ax.set_xlabel("Antioxidants")
    
