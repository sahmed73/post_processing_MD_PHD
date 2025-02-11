# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Feb 25 23:59:12 2024
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import os
import random
import sys
import pandas as pd
import magnolia.plot_template as mplt
from scipy.optimize import curve_fit

baseoil = 'PAO4'
location = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_6A_Soria\Production\1600'
sim_dir = ['Sim-1','Sim-2','Sim-3']    

atomsymbols= ['H','C','O']
bonddata = {}
for sim in sim_dir:
    bondfilepath = location+'\\'+sim+'\\bonds.reaxc'
    bonddata[sim]= bfp.parsebondfile(bondfilepath,cutoff=0.3)
#%%
timestep        = 0.25
ramp_rate       = 4
initial_temp    = 300
skipts          = (1600-initial_temp)/ramp_rate

df = pd.DataFrame()
for sim in sim_dir:
    neighbours = bonddata[sim]['neighbours']
    atypes    = bonddata[sim]['atypes']
    
    H_abs = {}
    for step, neigh in neighbours.items():
        molecules = bfp.get_molecules(neigh)
        ps = step*timestep/1000
        H_abs[ps] = 0
        for molecule in molecules:
            atom_types = [atypes[u] for u in molecule]
            nH = atom_types.count(1)
            nC = atom_types.count(2)
            if nC==30:
                H_abs[ps] += (62-nH)
                
    adf = pd.Series(H_abs)
    adf.index = adf.index-skipts
    adf = adf.loc[0:]
    
    df[sim]=adf
#%%
saveplt = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures'+f'\\H_abs_{random.randint(0,100000000)}'

fig, ax = plt.subplots()
# ax.set_plot(df,color='tab:red')
mplt.mean_range_plot(df, color='tab:red',ax=ax)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Number of H abstraction')
ax.set_ylim(top=40)
ax.set_title(location[location.find('Base'):]+'\n')
fig.savefig(saveplt,dpi=300,bbox_inches='tight')