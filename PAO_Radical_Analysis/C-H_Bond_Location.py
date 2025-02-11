# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Mar 14 02:15:46 2024
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

dirr      = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600"
sim_dirr  = ['Sim-1', 'Sim-2', 'Sim-3']

bonddata = {}
bo_cutoff = 0.3
for sim in sim_dirr:
    bondfilepath = dirr+'\\'+sim+'\\bonds.reaxc'
    bonddata[sim]=bfp.parsebondfile(bondfilepath, cutoff=bo_cutoff,
                                    bo=True)
#%%
ref = None
isomer_bonds_sim = {}
for sim in sim_dirr:
    bondfilepath = dirr+'\\'+sim+'\\bonds.reaxc'
    isomer_bonds, ref = bfp.map_isomer_bonds(bondfilepath,ref=ref,
                                             molecule_length=92)
    isomer_bonds_sim[sim]=isomer_bonds
a,b,c = isomer_bonds_sim.values()
isomer_bonds_flat = a+b+c
#%%
start = time.time()
first_CH = []
for sim in sim_dirr:
    neighbours = bonddata[sim]['neighbours']
    atypes     = bonddata[sim]['atypes']
    bondorders = bonddata[sim]['bondorders']
    
    isomer_bonds = isomer_bonds_sim[sim]
    for molecule_bonds in isomer_bonds:
        break_steps = []
        for i, (u, v) in enumerate(molecule_bonds):
            break_flag = False
            for step in bondorders:
                bo = bondorders[step][u].get(v,0)
                if bo<0.3:
                    break_steps.append(step)
                    break_flag = True
                    break
            if not break_flag:
                break_steps.append(step)
        
        min_step = np.inf
        for i in range(len(molecule_bonds)):
            u, v  = molecule_bonds[i]
            types = (atypes[u],atypes[v])
            if break_steps[i]<min_step and types in [(1,2),(2,1)]:
                min_step = break_steps[i]
                min_bond = molecule_bonds[i]        
        first_CH.append(min_bond)
        
end = time.time()
print(end-start)
#%%
# removing C-C bonds
CH_isomer_bonds_flat = [[] for i in range(75)]
for i, each in enumerate(isomer_bonds_flat):
    count = 0
    for u, v in each:
        if (atypes[u],atypes[v]) in [(1,2),(2,1)]:
            count+=1
            CH_isomer_bonds_flat[i].append((u,v))
            
bond_index = []
for f, iso in zip(first_CH,CH_isomer_bonds_flat):
    bond_index.append(iso.index(f))

arr = np.array(bond_index)
unique_elements, counts = np.unique(arr, return_counts=True)
# Create a dictionary from unique_elements and counts for better readability
element_counts = dict(zip(unique_elements, counts))

maxx = max(element_counts.values())
indices = []
for key, value in element_counts.items():
    if value==maxx:
        print(key)
        indices.append(key)
        
arrf = np.array(first_CH)
arrf[indices]






