# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 06:33:14 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1'
bondfilepath = directory+'\\bonds.reaxc'

bondinfo = bfp.parsebondfile(bondfilepath,bo=True)
#%%----
neighbours = bondinfo['neighbours']
atypes     = bondinfo['atypes']
bondorders = bondinfo['bondorders']
atomsymbols= ['H','C','O']

#%%
for step, neigh in neighbours.items():
    for parent, children in neigh.items():
        ptype  = atypes[parent]
        ctypes = [atypes[x] for x in children]
        if ptype==2 and len(ctypes)==3 and ctypes.count(3)==2 and ctypes.count(2)==1:
            print(parent,children)
            for child in children:
                bo=bondorders[step][parent][child]
                print("({},{})={}".format(atomsymbols[atypes[parent]-1],
                                          atomsymbols[atypes[child]-1],bo))
            print('----------------------------------')