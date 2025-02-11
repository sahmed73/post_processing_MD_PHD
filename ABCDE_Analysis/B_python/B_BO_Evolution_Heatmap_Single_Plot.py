# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 03:45:33 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne
import pandas as pd
import numpy as np
import seaborn as sns
import random

# directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\200_Less_TempRamp_1150K'
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\Pure_50B\Production\1050\Sim-1'
filename  = 'bonds.reaxc'
bondfilepath = directory+'\\'+filename

bondinfo = bfp.parsebondfile(bondfilepath,bo=True)
#%%------
import sys
neighbours = bondinfo['neighbours']
bondorders = bondinfo['bondorders']
atypes     = bondinfo['atypes']
steps      = list(neighbours.keys())
firststepneigh = neighbours[steps[0]]

bonds = {}
for parent,children in firststepneigh.items():
    #O-H bond
    key = 'O-H'
    if atypes[parent]==3:
        for child in children:
            if atypes[child]==1:
                bond = (parent,child)
                if key in bonds:
                    bonds[key].append(bond)
                else:
                    bonds[key]=[]
    
    #tert-butyl
    key = 'tert-butyl'
    if atypes[parent]==2:
        ff = [atypes[x] for x in children].count(2)
        if ff==4:
            for child in children:
                fff = [atypes[x] for x in firststepneigh[child]].count(2)
                if fff==3:
                    bond = (parent,child)
                    if key in bonds:
                        bonds[key].append(bond)
                    else:
                        bonds[key]=[]
    
    #link
    key = 'link'
    if atypes[parent]==2:
        ff = [atypes[x] for x in children].count(2)
        if ff==3:
            for child in children:
                fff = [atypes[x] for x in firststepneigh[child]].count(2)
                if fff==3:
                    bond = (parent,child)
                    if key in bonds:
                        bonds[key].append(bond)
                    else:
                        bonds[key]=[]
    #C-H
    key = 'C-H'
    if atypes[parent]==1:
        for child in children:
            if atypes[child]==2:
                bond = (parent,child)
                if key in bonds:
                    bonds[key].append(bond)
                else:
                    bonds[key]=[]
                    
    #C=C
    key = 'C=C'
    if atypes[parent]==2:
        Noxy = [atypes[x] for x in children].count(3)
        Nhyd = [atypes[x] for x in children].count(1)
        if Noxy==1 and Nhyd==2:
            for child in children:
                if atypes[child]==2:
                    for schild in firststepneigh[child]:
                        if atypes[schild]==2 and schild!=parent:
                            bond = (schild,child)
                            if key in bonds:
                                bonds[key].append(bond)
                            else:
                                bonds[key]=[]
    
    # -OH
    key = '-OH'
    if atypes[parent]==3:
        for child in children:
            if atypes[child]==2:
                bond = (parent,child)
                if key in bonds:
                    bonds[key].append(bond)
                else:
                    bonds[key]=[]

for b in bonds:
    print(b,len(bonds[b]))
#%%    
import magnolia.plot_template as mplt

plot_keys = {'O-H':'Bond-1 dissociation',
            '-OH':'Bond-2 dissociation',
            'link':'Bond-3 dissociation',
            'C=C':'Bond-4 dissociation',
            }

skipts = (1150-300)/4

fig, ax = plt.subplots(1,len(plot_keys)+1,figsize=(9, 4),
                       gridspec_kw={'width_ratios': [1, 1, 1, 1, 0.1]})
subtitle = ["O–H bond","C–O bond","C–C bond","C=C bond"]
for i,key in enumerate(plot_keys):
    b = bonds[key]
    steps = np.array(list(bondorders.keys()))
    timestep = 0.25
    df = bfp.bondorder_evolution(bondorders, b, ts=timestep,skipts=skipts,plot='no')

    hm = sns.heatmap(df,cmap='jet',xticklabels=1000,
                     ax=ax[i],vmin=0.0,vmax=2.0,cbar=False)
    
    xticklabels = [int(float(x.get_text())) for x in ax[i].get_xticklabels()]
    ax[i].set_xticklabels(xticklabels,rotation=90,fontsize=8)
    ax[i].set_yticks([])
    ax[i].set_title('{}B: {}'.format(i+1,subtitle[i]),fontsize=12)

cbar_ticks = [0,1,2]
cbar_kws = { 'label': 'Bond order','ticks': cbar_ticks}
sns.heatmap(df,cmap='jet',cbar_kws=cbar_kws,
            xticklabels=1000,ax=ax[-2],vmin=0.0,vmax=2.0,
            cbar_ax=ax[-1])

xticklabels = [int(float(x.get_text())) for x in ax[-2].get_xticklabels()]
ax[-2].set_xticklabels(xticklabels,rotation=90)
ax[-2].set_yticks([])
fig.text(0.5, -0.04, 'Time (ps)', ha='center', va='center')

title = directory[directory.find('ABCDE'):]+"Molecule D\n"
fig.suptitle(title, y=1.04, fontsize=10)
plt.grid(True)
mplt.saveplot(folder='BO_Evolution',name='bo_evolution')