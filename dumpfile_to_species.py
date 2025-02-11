# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Mar  3 04:11:00 2024
"""

import pandas as pd
import matplotlib.pyplot as plt
import random
import magnolia.speciesfile_parser as sfp
import magnolia.bondfile_parser as bfp
import magnolia.dumpfile_parser as dfp
import numpy as np
               
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\20_PAO4_80_O2_6A_Soria\Onset\Sim-1'

dumpfile = directory+'\\onset.lammpstrj'
dumpdata = dfp.parsedumpfile(dumpfile)
#%%
positions = {}
atomtypes = {}
for step, values in dumpdata.items():
    positions[step]={}
    for atom, atom_info in values.items():
        atomtype, mol, *pos = atom_info
        xyz = np.array([float(u) for u in pos])
        positions[step][atom] = xyz
        atomtypes[atom] = atomtype
#%%
cutoff = 2.21
dist = []
neigh = {}
for step in positions:
    for u in positions[step]:
        neigh[u]=[]
        for v in positions[step]:
            duv = positions[step][u]-positions[step][v]
            d   = np.linalg.norm(duv)
            if d<=cutoff:
                neigh[u].append(v)
            dist.append(d)
    break
molecules = bfp.get_molecules(neigh)