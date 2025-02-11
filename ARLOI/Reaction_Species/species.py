# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Nov 19 13:43:11 2024
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from collections import defaultdict

dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0001\Production\300-500K_TRR=1Kpps\Sim-1'
bondfile=dirr+r'\bonds.out'
atominfo=bfp.parsebondfile(bondfile,mtypes=True)
#%%
neighbours=atominfo['neighbours']
atypes=atominfo['atypes']
mtypes=atominfo['mtypes']

timestep=0.25
prev_molecules={}
reaction_lineage={}
for step, neigh in neighbours.items():
    time=step*timestep/1000
    molecules=list(bfp.get_molecules(neigh))
    for current_molecule in molecules:
        parent_molecules = []
        for prev_molecule in prev_molecules:
            if current_molecule & prev_molecule:
                parent_molecules.append(prev_molecule)
    prev_molecules=molecules.copy()
    print(parent_molecules)