# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:21:13 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import magnolia.dumpfile_parser_v1 as dfp
import matplotlib.pyplot as plt
import sys
import pickle
import os
import time
import random
import seaborn as sns
import numpy as np
import pandas as pd

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Onset\Sim-1'
filename    = 'bonds.reaxc'
bondfile_path    = directory+'\\'+filename

neighbours,atomtypes,bondorders,mtypes = bfp.get_neighbours(bondfile_path, bo=True, mtypes=True)
steps = list(neighbours.keys())
#%%--

bonds = []
atoms = []

firststepneigh = neighbours[steps[0]]

for parent,children in firststepneigh.items():
    if atomtypes[parent]==3 and len(children)==1 and atomtypes[children[0]]==2:
        bonds.append((parent,*children))
        atoms.append(parent)
        atoms.append(*children)
#%%
dumpfile = directory+'\\'+'onset.lammpstrj'

data = dfp.parsedumpfile(dumpfile)
#%%
print(len(data))

a1,a2 = bonds[0]
print(a1,a2)

steps,d = [],[]

for step,info in data.items():
    a1_pos = info[a1][1]
    a2_pos = info[a2][1]
    
    d.append(dfp.distance(a1_pos, a2_pos))
    steps.append(step)
timestep = 0.25        
ps = bfp.step2picosecond(steps, timestep)

plt.style.use('classic')
plt.xlabel('Time (ps)')
plt.ylabel('Distance (A)')
plt.scatter(ps,d,s=50,alpha=0.5, edgecolor='black', linewidth=1)
plt.savefig('python_outputs\\distance_plot', dpi=400,bbox_inches='tight')