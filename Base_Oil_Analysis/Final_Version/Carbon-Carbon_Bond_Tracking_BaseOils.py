# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Nov 15 22:35:40 2023
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1"
bondfilepath = directory+"\\bonds.reaxc"
bonddata = bfp.parsebondfile(bondfilepath,mtypes=True,bo=True)
neighbours = bonddata['neighbours']
atypes         = bonddata['atypes']
mtypes         = bonddata['mtypes']
bondorders     = bonddata['bondorders']
#%%
carbon_carbon_bonds = []
for step, neigh in neighbours.items():
    for parent, children in neigh.items():
        for child in children:
            ptype = atypes[parent]
            ctype = atypes[child]
            if ptype*ctype==4 and parent<child:
                carbon_carbon_bonds.append((parent,child))
    break
carbon_carbon_bonds.sort()
# bonds = []
# for i in range(0,len(carbon_carbon_bonds),29):
#     bonds.append(carbon_carbon_bonds[i])
base_oil = [x for x in range(1,3000) if mtypes[x]==1]