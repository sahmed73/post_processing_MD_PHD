# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 23 05:15:03 2023
"""

############################################
##
#   Heatmap of most abundance species
#   Abundance based on mean number of species
#   at the last 100(or anything you want) ps.
#   Basically, it show only the stable products not intermediates
#   Using species file
##
############################################

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import random
import pandas as pd
import magnolia.speciesfile_parser as sfp
import seaborn as sns
import numpy as np


# User Inputs
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\heatmap_of_species'

commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1000K"

# commmon = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1250'

### Very Important !!!!!!!!!!!!!!!!!!!!!!!!!!
# if isSpecefile==True: use species file, else use bond file 
isSpeciesFile = True

bo_cutoff = 0.45
flag = True
simulations = ['\\Sim-1']#,'\\Sim-2']#,'\\Sim-3']
for sim in simulations:
    directory = commmon+sim
    speciesfile_path = directory+'\\species.out' if bo_cutoff==0.3 else directory+f'\\species_{bo_cutoff}.out'
    bondfile_path    = directory+'\\bonds.reaxc'
    
    if isSpeciesFile:
        # using species file (too fast)
        current = sfp.get_species_count(speciesfile_path)
    else:
        # using bond file (too slow)
        current = bfp.get_species_count(bondfile_path,
                                        atomsymbols=['H','C','O'],
                                        cutoff=bo_cutoff)
    
    if flag:
        output_species = current
        flag = False
    else:
        output_species = output_species.add(current,fill_value=0)

## other wise plot doesn't end at 1000 ps
# if output_species.columns[0]!= 0:
#     print('Added 0 column')
#     output_species.insert(0, 0, output_species[output_species.columns[0]])
#%%
# User Inputs
timestep = 0.25
skipts   = 0#(1000-300)/4
nspecies = 10
tailps   = 400  #ps
lignin   = {"A":'H30C19O3',"B":'H34C26O4'}
exclude  = ['O2',lignin['A']]
n_xticks = 4  # how many xtick division in x axis i.e. 1+n_xticks ticklabels
nsim     = len(simulations)

tick_fontsize  = 14
label_fontsize = 18

inputfile_name= "species.out" if isSpeciesFile else "bonds.reaxc"
details  = ('\n\nThis heatmap is generated using {} file. Number of '
            'species are average\nover {} simulations.'
            'Species are sorted by mean species count of last {} ps.\n'
            'bo_cutoff {} is used here, '
            'First {} species are showing here.\n'.format(inputfile_name, nsim,
                                                          tailps, bo_cutoff,
                                                          nspecies))
title    = commmon[commmon.find('ABCDE'):]+ details



species = output_species.copy()/3 # copy is important
# species.replace(0,1e-20,inplace=True) #to avoid log(0) if use log scale
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
norm = None#LogNorm(vmin=1, vmax=species.max().max()) #log-scale or set None
fig, ax = plt.subplots(figsize=(8,5))
sns.heatmap(species,cmap='jet',xticklabels=xticks_freq, ax=ax,norm=norm,
            vmax=12)

cbar = ax.collections[0].colorbar
cbar.ax.set_ylabel('Number of species', fontsize=label_fontsize)
cbar_ticks = cbar.ax.get_yticklabels()
cbar.ax.set_yticklabels(cbar_ticks,fontsize=tick_fontsize)
# make the xticks integers
xticklabels = ax.get_xticklabels()
xticks_ = np.array([float(label.get_text()) for label in xticklabels])
xticks  = xticks_.astype(int)
yticks  = ax.get_yticklabels()
ax.set_xticklabels(xticks,fontsize=tick_fontsize)
ax.set_yticklabels(yticks,fontsize=tick_fontsize)
# ax.set_xlim(0,100*4)

ax.set_xlabel('Time (ps)',fontsize=label_fontsize)
plt.tick_params(left=False)
plt.title(title,fontsize=10)
plt.yticks(rotation=0)
plt.savefig(savedir+'\\hm_species_{}'.format(random.randint(0,9999999999)),dpi=1000, bbox_inches='tight')
plt.show()
