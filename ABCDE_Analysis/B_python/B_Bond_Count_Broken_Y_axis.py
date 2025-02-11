# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:21:43 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np

directory    = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\200_Less_TempRamp_1150K'
filename     = '\\bonds.reaxc'
bondfilepath = directory+filename
atomsymbols  = ['H','C','O']
bondinfo = bfp.parsebondfile(bondfilepath, bo=True)
#%%
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
#%%
import magnolia.needless_essential as ne
marker = ['s','^','x','o','>','*']  
color    = ['b','c','m','r','g','grey']  
cutoff = 1
tol    = 0.3
skipts = (1150-300)/4
markevery = 150

#Broke Axis Input
broken_from = 201
broken_to   = 1450

fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.15)

for i,key in enumerate(bonds.keys()):
    print(key,len(bonds[key]))
    # print(bonds[key])
    ps, nbonds = bfp.get_nbondsVStime(bondorders, atypes, bonds[key],
                                      skipts=skipts,cutoff=cutoff, tol=tol)
    
    ax1.plot(ps,nbonds,label=key,marker=marker[i],
                 markevery=markevery,markersize=4,linewidth=0.8,color=color[i])
    ax2.plot(ps,nbonds,label=key,marker=marker[i],
                 markevery=markevery,markersize=4,linewidth=0.8,color=color[i])

###-plot-####

ax1.set_ylim(0, broken_from)  # outliers only
ax2.set_yticks(range(broken_to,1600+1,50))
ax2.set_ylim(bottom=broken_to)  # most of the data

# hide the spines between ax and ax2
ax2.spines.bottom.set_visible(False)
ax1.spines.top.set_visible(False)
ax2.xaxis.tick_top()
ax2.tick_params(labeltop=False)  # don't put tick labels at the top
ax1.xaxis.tick_bottom()

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax2.plot([0, 1], [0, 0], transform=ax2.transAxes, **kwargs)
ax1.plot([0, 1], [1, 1], transform=ax1.transAxes, **kwargs)

fig.text(0.03, 0.5, 'Number of bonds', va='center', rotation='vertical',
         fontsize=12)

###-plot-####
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule B\n'
# plt.title(title)
#plt.title('Number of bond vs time {}'.format(bond_type))
ax2.legend(bbox_to_anchor=(1,1))
plt.xlabel('Time (ps)')
plt.savefig('..\\..\\python_outputs\\bondplot\\bondplot-'+ne.randstr(), dpi=300,bbox_inches='tight')