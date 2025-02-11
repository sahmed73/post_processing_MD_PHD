# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 23 05:15:03 2023
"""

############################################
##
#   Heatmap of most abundance species
#   Abundance based on mean number of species
#   at the last 100 ps.
#   Basically, it show only the stable products not intermediates
#   Using bonds file
##
############################################



import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import random
import pandas as pd
import magnolia.speciesfile_parser as sfp
import seaborn as sns
import numpy as np

savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\heatmap_of_species_speciesfile'

commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Squalane\25_Squalane_200_O2_Soria\Production\1600"

flag = True
for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    directory = commmon+sim
    bondfilepath = directory+'\\bonds.reaxc'
    current = bfp.get_species_count(bondfilepath,['H','C','O'])
    if flag:
        output_species = current
        flag = False
    else:
        output_species = output_species.add(current,fill_value=0)
#%%
# User Input
number_of_main = 25  
timestep = 0.25
cutoff   = 5
skipts   = (1600-300)/4
nspecies = 15
tailps   = 100  #ps
exclude  = ['O2','H62C30']
n_xticks = 4  # how many xtick division in x axis i.e. 1+n_xticks ticklabels 
title    = directory[directory.find('Base'):]+'  (using bonds.reaxc)\n\n'+ \
           'Species are sorted by mean species count of last {} ps. '\
               'Top {} species are showing here\n'.format(tailps,nspecies)



species = output_species.T.copy()   # copy is important
species.replace(0,1e-20,inplace=True)
species.columns = species.columns*timestep/1000-skipts
species    = species.loc[:,0:] # skip the ramping steps
species    = species.drop(exclude)

## mean of last {tailps} ps
means          = species.loc[:,species.columns.max()-tailps:].mean(axis=1)
sorted_index   = means.sort_values(ascending=False).index.tolist()
species        = species.loc[sorted_index]

species       = species.head(nspecies)
species.index = bfp.make_molecular_formula_latex(species.index,sort=True)

xticks_freq = int(species.columns.size/n_xticks)
norm = LogNorm(vmin=1, vmax=species.max().max())
fig, ax = plt.subplots(figsize=(8,0.3*nspecies))
sns.heatmap(species,cmap='jet',xticklabels=xticks_freq,
                 norm=norm, ax=ax)

# make the xticks integers
xticklabels = ax.get_xticklabels()
xticks = np.array([float(label.get_text()) for label in xticklabels])
ax.set_xticklabels(xticks.astype(int))

ax.set_xlabel('Time (ps)')
plt.tick_params(left=False)
plt.title(title,fontsize=7)
plt.yticks(rotation=0)
plt.savefig(savedir+'\\hm_species_{}'.format(random.randint(0,9999999999)),dpi=500, bbox_inches='tight')
