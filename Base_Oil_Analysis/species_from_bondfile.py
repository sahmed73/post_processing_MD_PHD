# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 23 05:15:03 2023
"""

import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
import random
import pandas as pd

commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Squalane\25_Squalane_200_O2_Soria\Production\1600"

flag = True
for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    directory = commmon+sim
    bondfilepath = directory+'\\bonds.reaxc'
    current = bfp.get_species_count(bondfilepath,['H','C','O'])
    if flag:
        summed_species = current
        flag = False
    else:
        summed_species = summed_species.add(current,fill_value=0)
#%%
number_of_main = 25  
cutoff = 5
skipts = (1600-300)/4

avg_species = summed_species/3

last_count = avg_species.iloc[-1,:].sort_values(ascending=False)
sorted_species = avg_species[:,last_count.index]
cutoff_species = sorted_species.loc[:,(sorted_species>=cutoff).any()]
drop_species = cutoff_species.drop(['H62C30','O2'],axis=1)
drop_species.index = drop_species.index*0.25/1000-skipts
skip_species = drop_species.loc[skipts:,:]
#%%

species = skip_species
time    = species.index

plt.style.use("classic")
plt.figure(figsize=(6.4, 4.8))
for column in species.columns:
    label = bfp.make_molecular_formula_latex([column],sort=True)
    plt.plot(time,species[column],label=label,linewidth=0.8)
    
plt.legend(fontsize=15,loc='upper left')
plt.xlabel('Time (ps)',fontsize=15)
plt.ylabel('Number of molecules',fontsize=15)
plt.title(directory[directory.find('Base'):]+'\n',fontsize=8)
plt.grid('on')
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\Species_from_bondfile'+'\\species_{}'.format(random.randint(0,10000000000)),dpi=400,bbox_inches='tight')
