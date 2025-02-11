# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Jan 14 04:21:15 2024
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import numpy as np
import seaborn as sns
import random
import magnolia.speciesfile_parser as sfp

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1250\Sim-1'
filename  = 'bonds.reaxc'
bondfilepath = directory+'\\'+filename

bo_cutoff = 0.50
bonddata = bfp.parsebondfile(bondfilepath,cutoff=bo_cutoff,
                             mtypes=True,bo=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
asyms = 'HCO'

bondlist = ['left_ring_O-H', 'right_ring_O-H','left_chain_O-H','right_chain_O-H',
            'left_ring_O-C','right_ring_O-C','left_chain_O-C','right_chain_O-C',
            'OH-H']


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
    elif value == 99:
        pass


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
markers = ['s','^']
cyan = (0, 0.85, 0.85)
colors    = ['b','b','r','r','m','m',cyan,cyan,'g']
cutoff = 1
tol    = 0.9
skipts = 0#(1150-300)/4
markevery = 150
tick_fontsize  = 14
label_fontsize = 18
markersize = 6
linewidth = 1.2
markeredgewidth = 1

species = sfp.get_species_count(directory+'\\species.out')
# species.index = species.index*0.25/1000
# species = species.loc[skipts:,:]
# species.index = species.index-skipts
print(species.loc['H2O'].values.size)
bonds['OH-H']=species.loc['H2O'].values[-4000:]/3
# print(bonds['OH-H'])

bond_number = ['1B','1B\'','2B','2B\'','3B','3B\'','4B','4B\'', 'HO-H bond formation']

fig, ax1 = plt.subplots(figsize=[8,6])
for i,key in enumerate(bondlist):
    if key in ['OH-H']:
        ax2 = ax1.twinx()
        ax2.plot(ps,bonds[key],label='HO—H',linewidth=linewidth,
                 color=colors[i], marker=markers[i%len(markers)],
                 markevery=214 if i==5 else 250, markersize=markersize,
                 markeredgewidth=markeredgewidth, markeredgecolor='black')
        ax2.set_ylim(0,52)
        ax2.set_ylabel("Number of bonds formed", fontsize=label_fontsize)
        ax2.legend(bbox_to_anchor=(1.56,0.68),title='Bond Formation',
                   fontsize=tick_fontsize, title_fontsize=label_fontsize)
        continue
    print(key,len(bonds[key]))
    ps, nbonds = bfp.get_nbondsVStime(bondorders, atypes, bonds[key],
                                      skipts=skipts,cutoff=cutoff, tol=tol)

    # plt.plot(ps,nbonds,label=key,marker=marker[i%len(marker)],
    #          markevery=markevery,markersize=4,linewidth=0.8,
    #          color=color[i%len(color)])
    ax1.plot(ps,nbonds,label=bond_number[i]+": "+key[-3:].replace('-','—'),
             linewidth=linewidth, color=colors[i],
             marker=markers[i%len(markers)],markevery=217 if i==5 else 250,
             markersize=markersize,
             markeredgewidth=markeredgewidth, markeredgecolor='black')

###-plot-####
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule B\n'
plt.title(title)
#plt.title('Number of bond vs time {}'.format(bond_type))
# Assume handles, labels are obtained from the plot
handles, labels = ax1.get_legend_handles_labels()

# Reorder handles and labels to match desired pattern
new_labels = [labels[i] for i in range(0, len(labels), 2)] + [labels[i] for i in range(1, len(labels), 2)]
new_handles = [handles[i] for i in range(0, len(handles), 2)] + [handles[i] for i in range(1, len(handles), 2)]
ax1.legend(new_handles, new_labels, bbox_to_anchor=(1.09,1.03),ncol=2,
           title= 'Bond Dissociation', fontsize=tick_fontsize,
           title_fontsize=label_fontsize)
ax1.set_xlabel('Time (ps)',fontsize=label_fontsize)
ax1.set_ylabel('Number of bonds dissociated',fontsize=label_fontsize)
ax1.tick_params(axis='both', which='major', labelsize=tick_fontsize)
ax2.tick_params(axis='both', which='major', labelsize=tick_fontsize)
plt.xlim(-1,1001)
# plt.ylim(-5,205)
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\bondplot'
# plt.grid('on')
plt.savefig(savedir+'//_B-'+str(random.randint(0, 999999999)), dpi=500,bbox_inches='tight')