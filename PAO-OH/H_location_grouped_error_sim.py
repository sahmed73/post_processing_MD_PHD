# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Sep 13 12:55:54 2024
"""

# H location

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne
import pandas as pd
import networkx as nx
import numpy as np
import seaborn as sns

data = {'Antioxidant': [], 'H location': [], 'Count': [], 'Error': []}

AOs = ['A0001', 'A0002', 'A0003', 'A0004', 'A0005']
for AO in AOs:
    # if AO == 'A0001':
    #     directory=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\1_ARLOI_24-08-28--00-09-01\15_PAO-OH_20_A0001\Production\Sim-2_with_H_removed'
    # else:
    #     directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_{}\Production\Sim-1'.format(AO)
    
    pos_count = []
    neg_count = []
    for sim in ['Sim-1','Sim-2','Sim-3']:
        directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_{}\Production\{}'.format(AO,sim)
        filename  = '\\bonds.reaxc'
        bondfile  = directory+filename
        speciesfile = directory+'species.out'
        atomsymbols = ['H','C','O']
        atominfo = bfp.parsebondfile(bondfile,bo=True)
        
        neighbours = atominfo['neighbours']
        atypes = atominfo['atypes']
        bondorders = atominfo['bondorders']
        asyms = 'HCO'
        
        H_location = []
        
        checked = []
        for step, neigh in neighbours.items():
            molecules = bfp.get_molecules(neigh)
            graph    = nx.Graph(neigh)
            for molecule in molecules:
                species = bfp.get_molecular_formula(molecule, atypes, asyms)
                if species=='H62C30O' and molecule not in checked:
                    checked.append(molecule)
                    subgraph = graph.subgraph(molecule).copy()
                    
                    for node in subgraph.nodes:
                        if atypes[node]==3:
                            connected_nodes = list(subgraph.neighbors(node))
                            for cnode in connected_nodes:
                                if atypes[cnode]==1:
                                    H_location.append(cnode)
                        
             
        positive=[]
        negative=[]
        
        for loc in H_location:
            for step, neigh in neighbours.items():
                graph = nx.Graph(neigh)
                connected_nodes = list(graph.neighbors(loc))
                for cnode in connected_nodes:
                    if atypes[cnode]==3:
                       positive.append(loc)
                    else:
                        negative.append(loc)
                break
        
        pos_count.append(len(positive))
        neg_count.append(len(negative))              
    
    mean_pos_count = np.mean(pos_count)
    mean_neg_count = np.mean(neg_count)
    
    data['Count'].append(mean_pos_count)
    data['Count'].append(mean_neg_count)
    data['Error'].append(np.std(pos_count))
    data['Error'].append(np.std(neg_count))
    data['H location'].extend(['Hydroxyl-H', 'Nonhydroxyl-H'])
    data['Antioxidant'].extend([AO,AO])
#%%
# Create a DataFrame
df = pd.DataFrame(data)

# Create the grouped bar plot
fig, ax = plt.subplots(dpi=300, figsize=[10,6])
# sns.barplot(x='Antioxidant', y='Count', hue='H location', yerr=data['Error'].values,
#             data=df, ax=ax)

sns.barplot(x='Antioxidant', y='Count', hue='H location', data=df, ci=None, ax=ax, capsize=0.1)
plt.xticks(rotation=45)

# Iterate over the container of bar objects
yerr=data['Error']
for i, bar in enumerate(ax.patches):  # ax.patches contains each bar
    ax.errorbar(
        bar.get_x() + bar.get_width() / 2,  # Get the x position of the bar
        bar.get_height(),  # Get the height (y value) of the bar
        yerr=yerr[i % len(yerr)],  # Error value for each bar
        fmt='none',  # No line connecting error bars
        capsize=5,  # Add caps to error bars
        color='black'  # Color of error bars
    )

# plt.subplots(dpi=350)
# plt.rcParams['font.size'] = 16
# plt.bar(x,y,width=0.5,color='maroon')
# plt.ylabel('Count')
