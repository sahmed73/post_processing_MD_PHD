# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Oct 27 11:19:49 2023
"""
############################################
##
#   Heatmap of most abundance species
#   Abundance based on mean number of species
#   at the last 100 ps.
#   Basically, it show only the stable products not intermediates
#   Using species file
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

baseoil = 'PAO4'

savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\heatmap_of_species'

commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600".format(baseoil,baseoil)

flag = True
cutoff = 0.30
for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    directory = commmon+sim
    speciesfile_path = directory+f'\\species_{cutoff}.out'
    if cutoff == 0.30: speciesfile_path = directory+'\\species.out'
    current = sfp.get_species_count(speciesfile_path)
    if flag:
        output_species = current
        flag = False
    else:
        output_species = output_species.add(current,fill_value=0)

#%%
# User Input
number_of_main = 25  
timestep = 0.25
skipts   = (1600-300)/4
nspecies = 15
tailps   = 10  #ps
exclude  = ['O2','H62C30']
n_xticks = 4  # how many xtick division in x axis i.e. 1+n_xticks ticklabels 
details  = ('\n\nThis heatmap is generated using species.out file.\nNumber of '
            'species are average over three simulations.\n'
            'Species are sorted by mean species count of last {} ps.\n'
            'First {} species are showing here.\n'.format(tailps, nspecies))
title    = commmon[commmon.find('Base'):]+ details+f"(cutoff: {cutoff})\n"

## plot parameter
label_fontsize  = 18
tick_fontsize   = 14



species = output_species.copy()/3 # copy is important
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
norm = None#LogNorm(vmin=1, vmax=60)

plt.style.use('default')
fig, ax = plt.subplots(figsize=[8,6])

heatmap = sns.heatmap(species,cmap='jet',xticklabels=xticks_freq,
                 ax=ax,norm=norm,vmax=50)

# make the xticks integers
xticklabels = ax.get_xticklabels()
xticks = np.array([float(label.get_text()) for label in xticklabels])
ax.set_xticklabels(xticks.astype(int))
ax.tick_params(axis='both', labelsize=tick_fontsize)

## Color Bar
cbar = heatmap.collections[0].colorbar  # Get the colorbar object
cbar.set_label('Number of species', fontsize=label_fontsize)  # Set the label and fontsize
# Customize other colorbar properties if needed
cbar.ax.tick_params(labelsize=tick_fontsize)  # Set tick label size

ax.set_xlabel('Time (ps)',fontsize=label_fontsize)
plt.tick_params(left=False)
plt.title(title,fontsize=6)
plt.yticks(rotation=0)
plt.savefig(savedir+'\\hm_species_{}'.format(random.randint(0,9999999999)),dpi=500, bbox_inches='tight')
