# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Mar  5 16:18:26 2024
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1'
bondfile_path = dirr+'\\bonds.reaxc'
bonddata = bfp.parsebondfile(bondfile_path,bo=True,mtypes=True)
#%%
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
mtypes     = bonddata['mtypes']
bondorders = bonddata['bondorders']
timestep   = 0.25
atomsymbols= 'HCO'

## only for molecule 1
molecule_id = 1
bo_evo = {}
for step, neigh in neighbours.items():
    ps = step*timestep/1000
    bo_evo[ps] = {}
    for parent, children in neigh.items():
        for child in children:
            if mtypes[parent]==molecule_id and mtypes[child]==molecule_id:
                if parent<child: continue
                key = f'{parent}-{child}'
                bo_evo[ps][key]=bondorders[step][parent][child]
#%%--plotting
df = pd.DataFrame(bo_evo).fillna(0).T
pdf = df.iloc[:,:]
plt.plot(pdf,label=pdf.columns)
# plt.legend(fontsize=10, loc='lower left')