# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Feb  6 04:44:19 2024
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne
import pandas as pd
import numpy as np
import seaborn as sns
import random

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1100K'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
atominfo = bfp.get_neighbours(bondfile,bo=True)
#%%--
neighbours, atypes, bondorders = atominfo
steps = list(neighbours.keys())
bondlist = []

firststepneigh = neighbours[steps[0]]
# 1 H
# 2 C
# 3 O

bonds = {'C-H':[],'O-H':[]}
for parent,children in firststepneigh.items():    
    if atypes[parent] == 1 and len(children)==1:
        child = children[0]
        if atypes[child]==2:
            bonds['C-H'].append((parent,child))
        elif atypes[child]==3:
            bonds['O-H'].append((parent,child))
        else:
            print('H-H')
            
broke = {'C-H':0,'O-H':0}

for key, value in bonds.items():
    for u, v in value:
        for step in bondorders.keys():
            bo = bondorders[step][u].get(v,0)
            if bo<0.3:
                if atypes[v]==2:
                    broke['C-H']+=1
                    break
                elif atypes[v]==3:
                    broke['O-H']+=1
                    break
#%%
print(broke)   
for key in broke.keys():
    print(key,len(bonds[key]))
    print(key, broke[key]*100/len(bonds[key]))    
#%%-
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution'

n_plot = len(bonds)
timestep = 0.25
skipts = (1200-300)/4
tick_fontsize  = 11
label_fontsize = 15
vmin   = 0
vmax   = 2

fig, ax = plt.subplots(1,n_plot+1,figsize=(9, 4),
                       gridspec_kw={'width_ratios': [1]*n_plot+[0.1] })
subtitle = ["O–H bond","C–C bond","C–O bond","C–O bond"]
for i,key in enumerate(bonds):
    df = bfp.bondorder_evolution(bondorders, bonds[key], ts=timestep,
                                 skipts=skipts,plot='no')
    
    if i<n_plot-1:
        hm = sns.heatmap(df,cmap='jet',xticklabels=2000,
                         ax=ax[i],vmin=0.0,vmax=2.0,cbar=False)
    else:
        cbar_ticks = np.arange(vmin,vmax+1,step=1,dtype=int)
        cbar_kws = {'ticks': cbar_ticks}
        sns.heatmap(df,cmap='jet',cbar_kws=cbar_kws,
                    xticklabels=2000,ax=ax[i],vmin=vmin,vmax=vmax,
                    cbar_ax=ax[-1])
        ax[-1].set_ylabel("Bond order", fontsize=label_fontsize)
    
    xticklabels = [int(float(x.get_text())) for x in ax[i].get_xticklabels()]
    ax[i].set_xticklabels(xticklabels,rotation=90,fontsize=tick_fontsize)
    yticks = [1]+list(range(10,51,10))
    ax[i].set_yticks(yticks)
    ax[i].tick_params(bottom=True, left=True, right=False, top=False)
    if i==0: ax[i].set_yticklabels(yticks,rotation=0,fontsize=tick_fontsize)
    ax[i].set_title('{}A: {}'.format(i+1,subtitle[i]),
                    fontsize=label_fontsize-1)

fig.text(0.5, -0.04, 'Time (ps)', ha='center', va='center',
         fontsize=label_fontsize)
fig.text(0.08,0.5, 'Molecule Index', ha='center', va='center',rotation=90,
         fontsize=label_fontsize)


title = directory[directory.find('ABCDE'):]
fig.suptitle(title, y=1.04, fontsize=12)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution\bo_evolution_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')
