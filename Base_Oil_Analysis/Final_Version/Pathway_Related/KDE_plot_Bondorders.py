# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Feb  2 03:33:59 2024
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import numpy as np
import seaborn as sns
import random
from collections import Counter

baseoil = "Squalane"
simdir = ['Sim-1','Sim-2','Sim-3']

bonddata = {}
cutoff   = 0.35
for sim in simdir:
    directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(baseoil,baseoil,sim)
    bondfilepath = directory+'\\bonds.reaxc'
    bonddata[sim] = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%%
sim = 'Sim-1'

neighbors = bonddata[sim]['neighbours']
atypes    = bonddata[sim]['atypes']
mtypes    = bonddata[sim]['mtypes']
bondorders= bonddata[sim]['bondorders']
atomsymbols = 'HCO'


atom_type_pairs = [[1,2],[1,3],[2,2],[2,3],[3,3]]
atom_symb_pairs = ['C-H','O-H','C-C','C-O','O-O']
pairs = zip(atom_type_pairs,atom_symb_pairs)

data = {}
for check, symb in pairs:
    data[symb] = []
    print(check,symb)
    for step in bondorders:
        for u in bondorders[step]:
            for v in bondorders[step][u]:
                pair = sorted([atypes[u],atypes[v]])
                bo = bondorders[step][u][v]
                if check==pair:
                    data[symb].append(bo)

#%%
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\bond_order_KDE_{}'.format(random.randint(0,100000000))

label_fontsize = 20 
for key in data:
    print(key)   
    sns.kdeplot(data[key], label=key, fill=True)
    
plt.legend()
plt.xlabel("Bond order")
plt.minorticks_on()
plt.grid(True)
plt.savefig(savedir, dpi=500,bbox_inches='tight')