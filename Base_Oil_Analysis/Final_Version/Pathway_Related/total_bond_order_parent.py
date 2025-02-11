# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Feb  2 04:28:21 2024
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
import pandas as pd
import random

baseoil = "PAO4"
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

data = {}
for frame, (step, neigh) in enumerate(neighbors.items()):
    data[frame]={}
    for parent, children in neigh.items():
        total = 0
        for child in children:
            bo = bondorders[step][parent][child]
            total+=bo
        data[frame][parent]=total
#%%
atypes_series = pd.Series(atypes)
df  = pd.DataFrame(data)
df.columns = df.columns-1300
df = df.loc[:,df.columns>=0]
df = df.sort_index()
df.columns = df.columns*0.25

#  hydrogen
dfh = df[df.index.map(atypes) == 1]
dfh = dfh[dfh.gt(1.).any(axis=1)]

dfc = df[df.index.map(atypes) == 2]
dfo = df[df.index.map(atypes) == 3]
#%%
plt.figure(figsize=[10,5])
sns.heatmap(dfh, cmap='jet',vmax=2.0,vmin=0.0,
            cbar_kws={"label": "Total bond order"})
plt.xlabel("Time (ps)")
plt.ylabel("H atoms")
plt.yticks([])
plt.xticks(np.array([0, 250, 500, 750, 1000])*4, labels=['0', '250', '500', '750', '1000'],fontsize=8,rotation=0)
plt.savefig(r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution\total_bond_H_{}".format(random.randint(0,100000000)),dpi=500,bbox_inches='tight')
plt.show()