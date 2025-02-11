# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Nov 21 00:54:51 2023
"""

from collections import deque
import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man

directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1"
bondfilepath = directory+'\\bonds.reaxc'
bonddata = bfp.parsebondfile(bondfilepath,mtypes=True)
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
#%%

## getting all pao4 or squalane molecule's initials structure
number_of_baseoils = 25
molecule_adjLists = [{} for i in range(number_of_baseoils)]
for parent,children in neighbors[0].items():
    if atypes[parent]==1:
        continue
    
    if mtypes[parent] in range(1,number_of_baseoils+1):
        #oignoring the hydrogen
        children_without_H = [x for x in children if atypes[x]==2]
        molecule_adjLists[mtypes[parent]-1][parent]=children_without_H


frames = []
iii= 0 
for molecule_adjList in molecule_adjLists:
    iii+=1
    print(iii)
    frames.append(bfp.frame_when_bond_breaks(neighbors, molecule_adjList))
#%%
import numpy as np
frames_array = np.array(frames)
frames_array[frames_array==None] = 0
plt.scatter(range(25),sorted(frames_array),c='r')