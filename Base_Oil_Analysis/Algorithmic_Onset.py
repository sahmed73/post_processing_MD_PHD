# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 22:31:32 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\Onset\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
bondinfo = bfp.parsebondfile(bondfile,mtypes=True)
#%%--
neighbours = bondinfo['neighbours']
atypes     = bondinfo['atypes']
mtypes     = bondinfo['mtypes']
steps      = list(neighbours.keys())

number_of_molecule = 25

f = open('weiredABCDE.txt','w')

association = {}  # association[step][moldID]=[atomID_1, atomID_2,......]
for step, neigh in neighbours.items():
    association[step]={key:[] for key in range(1,number_of_molecule+1)}
    for parent, children in neigh.items():
        parentMolID = mtypes[parent] # molecule ID of parent
        if parentMolID in range(1,number_of_molecule+1):
            for child in children:
                childMolID = mtypes[child] # molecule ID of child
                if parentMolID != childMolID:
                    # (-) sign represent foreign atom association
                    if -child not in association[step][parentMolID]:
                        association[step][parentMolID].append(-child)
                else:
                    if child not in association[step][parentMolID]:
                        association[step][parentMolID].append(child)
count = {}
for step, ass in association.items():
    picosecond = 300+(step*0.25/1000)*4
    count[picosecond]=0
    print('Timestep:', step,'{}ps'.format(picosecond),file=f)
    for mols, atoms in ass.items():
        types = [atypes[abs(x)] for x in atoms]
        H_count = types.count(1)
        C_count = types.count(2)
        if C_count==30 and H_count==62 and len(types)==92:
            count[picosecond]+=1
        print(mols,'length',len(atoms),C_count,H_count,file=f)
        print(atoms,file=f)
        print(types,file=f)

plt.plot(count.keys(),count.values())
f.close()
