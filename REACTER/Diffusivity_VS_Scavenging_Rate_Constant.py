# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Mar 12 14:50:00 2025
"""

import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AOs          = ['A0001', 'A0002', 'A0003', 'A0004', 'A0005']
MAP = dict(zip(AOs, ['L1','L2','S1','S2','S3']))
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
    
#%% MSD plot
colors = ['tab:blue', 'tab:orange','tab:green','tab:red','tab:purple']
plt.style.use('default')
plt.rcParams['font.size']=18
fig, ax = plt.subplots(dpi=350)
for col in msd_data:
    if col=='Time': continue
    ax.plot(time/1000, msd_data[col], label=MAP[col], linestyle='--', linewidth=2)

ax.set_xlabel('Time (ns)')
ax.set_ylabel('MSD ($Å^2$)')
ax.legend()
#%% Diffusivity plot
dfs = []
for col in msd_data:
    if col=='Time': continue
    m, c = np.polyfit(msd_data['Time']/1000, msd_data[col], 1)
    D = m/6
    dfs.append(D)

# Plot bar chart for diffusivities
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
fig, ax = plt.subplots(dpi=350)
ax.bar([MAP[AO] for AO in AOs], dfs, width=0.7, color=colors, edgecolor='k', linewidth=1)
ax.set_ylabel("Diffusivity ($Å^2/ns$)")
ax.set_xlabel("Antioxidants")
ax.set_ylim(0,15.5)
ax.set_yticks(range(0,16,3))

#%% D vs K
##================
fig, ax2 = plt.subplots(dpi=350)
data = {
    "Sim-1": [0.081444, 0.077283, 0.074504, 0.080203, 0.116402],
    "Sim-2": [0.076825, 0.122914, 0.046360, 0.045120, 0.082771],
    "Sim-3": [0.044220, 0.098843, 0.058524, 0.077844, 0.096224]
}

df = pd.DataFrame(data, index=AOs)
k = df.mean(axis=1)

# Apply glowing effect using multiple scatter layers
for alpha, size in zip([0.1, 0.2, 0.3, 0.4, 0.5], [300, 250, 200, 150, 100]):
    ax2.scatter(k, dfs, s=size, color=colors, alpha=alpha,
                edgecolor='k')

ax2.set_xlabel("Scavenging Rate Constant (ns$^{-1}$)")
ax2.set_ylabel("Diffusivity ($Å^2/ns$)")
ax2.set_xlim(left=0.055, right=0.105)
ax2.set_ylim(0,15.5)
ax2.set_yticks(range(0,16,3))
ax2.grid(alpha=0.4, linestyle='--')

for i, txt in enumerate(AOs):
    txt = MAP[txt]
    if txt in []: offset = (-8,-18)
    else: offset = (5,5)
    ax2.annotate(txt, 
                 (k[i], dfs[i]),  # Position
                 fontsize=14, 
                 ha='left', va='bottom',  # Adjust alignment
                 xytext=offset,  # Offset text slightly (x, y)
                 textcoords='offset points')
    
