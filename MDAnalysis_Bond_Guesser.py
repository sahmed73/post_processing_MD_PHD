# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Oct 24 16:08:35 2024
"""

import MDAnalysis as mda
import numpy as np

path_to_data = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\9_ARLOI_Bond_OFF_High_Temp_Press\A0001\Equilibration\Sim-1"
u = mda.Universe(path_to_data + "\eq_1.data",
                 path_to_data + "\eq.npt-high.unwrapped.lammpstrj",
                 topology_format="data", format="lammpsdump",
                 atom_style='id mol type x y z',
                 guess_bonds=True, vdwradii={'1': 1.2, '2': 1.7, '3': 1.52})
