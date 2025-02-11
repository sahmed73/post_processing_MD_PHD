# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec  8 13:29:43 2023
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import magnolia.speciesfile_parser as sfp
import seaborn as sns
               
base_oil  = 'Squalane'
real_name = dict(PAO4='PAO 4', Squalane='Squalane')
n_distinct_species = {}
dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600".format(base_oil,base_oil)

sim_dir = ['Sim-1','Sim-2','Sim-3']
for sim in sim_dir:
    n_distinct_species[sim] = {}
    for cutoff in np.array(range(30,76,5))/100:
        directory = dirr+'\\'+sim
        speciesfile = directory+f'\\species_{cutoff}.out'
        if cutoff==0.30:
            speciesfile = directory+'\\species.out'
            
        species_count = sfp.get_species_count(speciesfile).T
        n_distinct_species[sim][cutoff] = species_count.columns.size
        
n_distinct_species = pd.DataFrame(n_distinct_species)
#%%
colors = ['r','b','g']
linewidth = 1.0
marker = 's'
markersize = 5
markerfacecolor = 'none'
markeredgecolor = colors
markeredgewidth = 1.5
label_fontsize  = 18
tick_fontsize   = 14

plt.style.use('default')
fig, ax = plt.subplots(figsize=[8,6])

x = n_distinct_species.index  
minn = n_distinct_species.min()
for i,column in enumerate(n_distinct_species.columns):
    y = n_distinct_species[column]
    ax.plot(x,y,label=column, marker=marker, markersize=markersize,
            markerfacecolor=markerfacecolor, color = colors[i],
            markeredgecolor=markeredgecolor[i], linewidth=linewidth,
            markeredgewidth=markeredgewidth)
    minx, miny = y.idxmin(),y.min()
    ax.plot(minx,miny,color=colors[i],marker='X',markersize=12,
            markeredgecolor=markeredgecolor[i])
    
ax.legend(loc='upper left',title = real_name[base_oil],
          title_fontsize=label_fontsize, fontsize = tick_fontsize)

ax.set_xlabel('Bond order cutoff',fontsize=label_fontsize)
ax.set_ylabel('Total number of distinct species',fontsize=label_fontsize)
ax.tick_params(axis='both', labelsize=tick_fontsize)
# ax.tick_params(axis='x', labelsize=12)
ax.set_xlim(0.29,0.8)
ax.set_ylim(np.array(ax.get_ylim())+np.array([-20,20]))
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+'\\species_speciesfile{}'.format(random.randint(0,10000000000)),dpi=500,bbox_inches='tight')