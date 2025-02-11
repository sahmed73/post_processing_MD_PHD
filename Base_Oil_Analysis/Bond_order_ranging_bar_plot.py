# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec  1 15:58:36 2023
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man

baseoil = "PAO4"
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(baseoil,baseoil)
bondfilepath = directory+'\\bonds.reaxc'
cutoff = 0.45
bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
atomsymbols = 'HCO'

boflat = []
for step, __ in bondorders.items():
    for atom1, _ in __.items():
        for atom2, bo in _.items():
            boflat.append(bo)
#%%
import numpy as np
data = np.array(boflat)
ranges = [0.45, 0.6, 1.1, 1.4, 2]

# Count elements in each range
counts = []
for i in range(len(ranges) - 1):
    lower_bound = ranges[i]
    upper_bound = ranges[i + 1]
    # Count how many data elements are in this range
    count_in_range = ((data >= lower_bound) & (data < upper_bound)).sum()
    counts.append(count_in_range)
    
counts_ = np.array(counts)/(1000)

# Bar plot
plt.bar(range(len(counts_)), counts_, tick_label=[1,2,3,4],color='tab:red')
plt.xlabel('Ranges')
plt.ylabel('Count')
plt.title('Bar Plot of Data in Different Ranges')
plt.show()
