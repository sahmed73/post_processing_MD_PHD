# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Apr  4 10:08:37 2024
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO_Radical\20_PAO_Radical_20_A_Soria\Production\800K\Sim-1'
filename = '\\bonds.reaxc'
bondfile = dirr+filename
bonddata = bfp.parsebondfile(bondfile,cutoff=0.3)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
steps     = np.array(list(neighbors.keys()))
ps = bfp.step2ps(steps,timestep=0.25)

# counts = []
# steps  = []
# for step, neigh in neighbors.items():
#     molecules = bfp.get_molecules(neigh)
#     count = 0
#     for molecule in molecules:
#         if len(molecule)==93:
#             count+=1
#     # print(count,end=',')
#     counts.append(count)
#     steps.append(step)
# plt.scatter(ps,counts)
# plt.show()
#%%
rad  = np.array([92+x*92 for x in range(20)])
H_OH = np.array([1865+x*52 for x in range(20)])
Habs_counts = []
for step, neigh in neighbors.items():
    count = 0
    for parent, children in neigh.items():
        if parent in rad:
            common = set(children) & set(H_OH)
            # print(len(common))
            if len(common)==1:
                count+=1
    Habs_counts.append(count)
    
with open('temporary.txt', 'a') as f:
    for number in Habs_counts:
        f.write(str(number) + "\n")
#%%
plt.plot(ps,Habs_counts)