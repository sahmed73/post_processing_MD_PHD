# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Jan  9 06:05:42 2024
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import numpy as np
import seaborn as sns
import random

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\Sim-1'
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

bondlist = ['left_ring_O-H', 'right_ring_O-H','left_chain_O-H','right_chain_O-H',
            'left_ring_O-C','right_ring_O-C','left_chain_O-C','right_chain_O-C',]
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
    
    left, right = sorted(center)
    
    for ox in oxygen:
        dist = nx.dijkstra_path_length(G, source=ox, target=left)
        length[ox]=dist
    

oxy_selection = {}
oxy_type_list = ['right_ring_O','left_ring_O','right_chain_O','left_chain_O']
for ox in oxy_type_list:
    oxy_selection[ox]=[]
    
for key, value in length.items():
    if value == 2:
        oxy_selection['left_ring_O'].append(key)
    elif value == 3:
        oxy_selection['right_ring_O'].append(key)
    elif value == 6:
        oxy_selection['left_chain_O'].append(key)
    elif value == 7:
        oxy_selection['right_chain_O'].append(key)


node_to_remove = [node for node in G.nodes if atypes[node]==1]
G.remove_nodes_from(node_to_remove)
c = []
for node in G.nodes:
    if node in oxy_selection['left_ring_O']:
        c.append('r')
    elif node in oxy_selection['right_ring_O']:
        c.append('pink')
    elif node in oxy_selection['left_chain_O']:
        c.append('green')
    elif node in oxy_selection['right_chain_O']:
        c.append('lightgreen')
    else:
        c.append('tab:blue')

    
nx.draw_spring(G,node_color=c,with_labels=False)
draw_one=False   
    
for parent, children in neighbors[0].items():
    for ox in oxy_type_list:
        if parent in oxy_selection[ox]:
            for child in children:
                if atypes[child]==1:
                    bonds[ox+'-H'].append((parent,child))
                elif atypes[child]==2:
                    bonds[ox+'-C'].append((parent,child))
#%%
skipts = 0#(1150-300)/4
vmin   = 0
vmax   = 2
fontsize = 7

count = 1
n_plot = len(bonds)
fig, ax = plt.subplots(1,n_plot+1,figsize=(1.25*n_plot, 4),
                       gridspec_kw={'width_ratios': [1]*n_plot+[0.1]})
j = 0
for i,key in enumerate(bonds):
    b = bonds[key]
    timestep = 0.25
    df = bfp.bondorder_evolution(bondorders, b, ts=timestep,
                                 skipts=skipts,plot='no')
    
    if j<n_plot-1:
        hm = sns.heatmap(df,cmap='jet',xticklabels=2000,
                         ax=ax[j],vmin=vmin,vmax=vmax,cbar=False)
    else:
        cbar_ticks = np.arange(vmin,vmax+1,step=1,dtype=int)
        cbar_kws = { 'label': 'Bond order','ticks': cbar_ticks}
        sns.heatmap(df,cmap='jet',cbar_kws=cbar_kws,
                    xticklabels=2000,ax=ax[j],vmin=vmin,vmax=vmax,
                    cbar_ax=ax[-1])
    
    xticklabels = [int(float(x.get_text())) for x in ax[j].get_xticklabels()]
    ax[j].set_xticklabels(xticklabels,rotation=90,fontsize=fontsize)
    yticks = [1]+list(range(10,51,10))
    ax[j].set_yticks(yticks)
    ax[j].tick_params(bottom=True, left=True, right=False, top=False)
    if j==0: ax[j].set_yticklabels(yticks,rotation=0,fontsize=fontsize)
    
    ## title
    if 'left_' in key:
        title = '{}B\n{}'.format(count,key.replace('_',' '))
    else:
        title = '{}B\'\n{}'.format(count,key.replace('_',' '))
        count+=1
        
    ax[j].set_title(title,fontsize=fontsize+1.3)
    j+=1


fig.text(0.5, -0.02,'Time (ps)',ha='center',va='center',fontsize=fontsize*2)
fig.text(0.08,0.5, 'Molecule Index', ha='center', va='center',rotation=90,
         fontsize=fontsize*2)

title = directory[directory.find('ABCDE'):]+"Molecule B_{}\n".format(cutoff)
fig.suptitle(title, y=1.03, fontsize=fontsize//2)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution\bo_evolution_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')