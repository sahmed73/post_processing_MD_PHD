# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec 13 01:21:49 2024
"""

import MDAnalysis as mda
import numpy as np

path = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0003"
topofile = path+r"\DataFile\Bulk\100_PAOr_50_A0003_density=0.2.data"
trajfile = path+r"\Reaction\Sim-13_T=1000K\reaction.unwrapped.lammpstrj"


u = mda.Universe(topofile, trajfile,
                 topology_format="data", format="lammpsdump",
                 atom_style='id resid type x y z')

#%%

for resid in range(1,151,1):
    # Access the first frame
    u.trajectory[0]  # Sets the universe to the first frame
    num_atoms_first = len(u.select_atoms(f"resid {resid}"))
    # Correct way to get the time of the current frame
    time_first = u.trajectory.time
    print(f"First Frame: Time (ps): {time_first}, Atoms in resid {resid}: {num_atoms_first}")
    
    # Access the last frame
    u.trajectory[-1]  # Sets the universe to the last frame
    num_atoms_last = len(u.select_atoms(f"resid {resid}"))
    # Correct way to get the time of the current frame
    time_last = u.trajectory.time
    print(f"Last Frame: Time (ps): {time_last}, Atoms in resid {resid}: {num_atoms_last}")
    print('--------------------------------------------------\n')