# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Dec 27 15:53:32 2023
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import numpy as np
import seaborn as sns
import random

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1150\Sim-1'
filename  = 'bonds.reaxc'
bondfilepath = directory+'\\'+filename

cutoff = 0.30
bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
asyms = 'HCO'

bondlist = ['ring_OH','chain_OH','ring_CO','chain_CO']
bonds = {}
for key in bondlist:
    bonds[key] = []

number_intial_molecules = 50
length = {}
draw_one = True
for m in range(1,number_intial_molecules+1):
    g = {}
    for parent, children in neighbors[0].items():
        if mtypes[parent]==m:
            g[parent]=children.copy()
    G = nx.Graph(g)
    center = nx.center(G)
    oxygen = [o for o in G.nodes if atypes[o]==3]
    
    for ox in oxygen:
        length[ox]=0
        for c in center:
            dist = nx.dijkstra_path_length(G, source=ox, target=c)
            length[ox]+=dist
    
    # bonds['C-C_link'].append(tuple(center))
    

ring_OH  = [key for key,value in length.items() if value==min(length.values())]
chain_OH = [key for key,value in length.items() if value==max(length.values())]

node_to_remove = [node for node in G.nodes if atypes[node]==1]
G.remove_nodes_from(node_to_remove)
c = []
for node in G.nodes:
    if node in ring_OH:
        c.append('r')
    elif node in chain_OH:
        c.append('pink')
    else:
        c.append('tab:blue')

    
nx.draw_spring(G,node_color=c,with_labels=False)
draw_one=False   
    
for parent, children in neighbors[0].items():
    if parent in ring_OH:
        for child in children:
            if atypes[child]==1:
                bonds['ring_OH'].append((parent,child))
            elif atypes[child]==2:
                bonds['ring_CO'].append((parent,child))
                
    elif parent in chain_OH:
        for child in children:
            if atypes[child]==1:
                bonds['chain_OH'].append((parent,child))
            elif atypes[child]==2:
                bonds['chain_CO'].append((parent,child))
#%%
skipts = 0#(1150-300)/4
vmin   = 0
vmax   = 2

fig, ax = plt.subplots(1,len(bonds)+1,figsize=(2.25*len(bonds), 4),
                       gridspec_kw={'width_ratios': [1]*len(bonds)+[0.1]})
for i,key in enumerate(bonds):
    b = bonds[key]
    timestep = 0.25
    df = bfp.bondorder_evolution(bondorders, b, ts=timestep,
                                 skipts=skipts,plot='no')

    hm = sns.heatmap(df,cmap='jet',xticklabels=1000,
                     ax=ax[i],vmin=vmin,vmax=vmax,cbar=False)
    
    xticklabels = [int(float(x.get_text())) for x in ax[i].get_xticklabels()]
    ax[i].set_xticklabels(xticklabels,rotation=90,fontsize=8)
    ax[i].set_yticks([])
    ax[i].set_title('{}B: {}'.format(i+1,key),fontsize=10)

cbar_ticks = np.arange(vmin,vmax+1,step=1,dtype=int)
cbar_kws = { 'label': 'Bond order','ticks': cbar_ticks}
sns.heatmap(df,cmap='jet',cbar_kws=cbar_kws,
            xticklabels=1000,ax=ax[-2],vmin=vmin,vmax=vmax,
            cbar_ax=ax[-1])

xticklabels = [int(float(x.get_text())) for x in ax[-2].get_xticklabels()]
ax[-2].set_xticklabels(xticklabels,rotation=90)
ax[-2].set_yticks([])
fig.text(0.5, -0.02, 'Time (ps)', ha='center', va='center')

title = directory[directory.find('ABCDE'):]+"Molecule B_{}\n".format(cutoff)
fig.suptitle(title, y=0.98, fontsize=5)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution\bo_evolution_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')