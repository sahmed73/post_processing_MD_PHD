# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Sep  6 15:37:25 2024
"""

import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import os

AOs = ['A0001', 'A0003']
tsd = {'A0001': '1.80', 'A0003': '1.95'}
sim = {'A0001': 7, 'A0003': 6}

msd_data = {}
diffusivity_data = {}

for AO in AOs:    
    # Load the topology and trajectory files
    # Make sure to specify the correct format for LAMMPS dump
    topo = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\{AO}\DataFile\Bulk\100_PAOr_50_{AO}_density=0.2.data"
    trej = rf'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\{AO}\Reaction\Sim-{sim[AO]}_TSD={tsd[AO]}'
    u = mda.Universe(topo, trej, format='LAMMPSDUMP')
    
    # Select all atoms or a specific group of atoms
    all_atoms = u.select_atoms('all')
    
    # Select all atoms (you can modify this to select specific atoms or groups)
    AO_atoms = all_atoms[all_atoms.resids > 100]
    
    # Compute MSD using the Einstein relation
    msd = EinsteinMSD(AO_atoms, msd_type='xyz', fft=True, unwarp=True).run()
    
    # Access the MSD timeseries (MSD values for each timestep)
    msd_timeseries = msd.results.timeseries
    
    # Get the corresponding time steps
    time_step_interval = 0.25  # Time interval per frame in picoseconds (adjust based on your simulation settings)
    timesteps = np.array([ts.frame * time_step_interval for ts in u.trajectory])
    
    # Store the timesteps and MSD timeseries in the dictionary
    msd_data[AO] = (timesteps, msd_timeseries)
    
    # Calculate diffusivity from the slope of MSD vs. time in the linear region
    slope, intercept, r_value, p_value, std_err = linregress(timesteps[(timesteps > 80) & (timesteps < 280)],
                                                             msd_timeseries[(timesteps > 80) & (timesteps < 280)])
    diffusivity = slope / (2 * 3)  # Divide by 2d (where d=3 for 3D diffusion)
    
    # Store the diffusivity
    diffusivity_data[AO] = diffusivity
    print(f"Diffusivity for {AO}: {diffusivity:.4e} Å^2/ps")

#%% Now you can plot the data independently
plt.rcParams['font.size']=16
fig, ax = plt.subplots(dpi=350)
for AO, (timesteps, msd_timeseries) in msd_data.items():
    ax.plot(timesteps[(timesteps>50) & (timesteps<200)],
            msd_timeseries[(timesteps>50) & (timesteps<200)], label=f"{AO}")
    
ax.set_xlabel("Time (ps)")
ax.set_ylabel("MSD ($Å^2$)")
ax.set_ylim(0)
# ax.title("Mean Squared Displacement vs Time")
plt.grid(False)
ax.legend()

#%% Plot diffusivity as a bar chart
fig, ax = plt.subplots(dpi=350)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

AOs = list(diffusivity_data.keys())
diffusivities = list(diffusivity_data.values())

ax.bar(AOs, diffusivities, color=colors)
ax.set_xlabel("Antioxidant")
ax.set_ylabel("Diffusivity ($Å^2/ps$)")
# ax.set_title("Diffusivity of Each Antioxidant")
plt.grid(False)
plt.show()
