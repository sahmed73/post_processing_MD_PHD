# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Mar  3 14:28:20 2024
"""

import MDAnalysis as mda
import sys
import math
import matplotlib.pyplot as plt
from MDAnalysis.lib.distances import distance_array


# Load the trajectory and topology files

trajectory_file = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\20_PAO4_80_O2_6A_Soria\Onset\Sim-1\onset.lammpstrj'

topology_file = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\20_PAO4_80_O2_6A_Soria\DataFile\20_PAO4_80_O2_6A.data'

# Create a Universe object
u = mda.Universe(topology_file, trajectory_file, format='LAMMPSDUMP', dt=0.25)


# Loop over all timesteps in the trajectory
for ts in u.trajectory:
    print(f"Time step: {ts.frame}")
    for atom in u.atoms:
        print(f"Atom ID: {atom.id}, Position: {atom.position}")
    
    # If you only need positions of a specific group of atoms, you can do:
    group = u.select_atoms('resid 1')  # Example: Selecting water molecules
    for atom in group:
        print(f"Atom ID: {atom.id}, Position: {atom.position}")
    
    break

    