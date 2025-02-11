# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Sep  2 11:08:10 2024
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

data = {'Antioxidant': [], 'H location': [], 'Count': []}

AOs = ['A0001', 'A0002']#, 'A0003', 'A0004', 'A0005']
for AO in AOs:
    if AO == 'A0001':
        directory=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\1_ARLOI_24-08-28--00-09-01\15_PAO-OH_20_A0001\Production\Sim-2_with_H_removed'
    else:
        directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_{}\Production\Sim-1'.format(AO)
    
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
        
    data['Count'].append(len(positive))
    data['Count'].append(len(negative))
    data['H location'].extend(['Hydroxyl-H', 'Nonhydroxyl-H'])
    data['Antioxidant'].extend([AO,AO])

#%%
# Create a DataFrame
df = pd.DataFrame(data)

# Create the grouped bar plot
fig, ax = plt.subplots(dpi=300, figsize=[10,6])
sns.barplot(x='Antioxidant', y='Count', hue='H location', data=df, ax=ax)
plt.xticks(rotation=45)

# plt.subplots(dpi=350)
# plt.rcParams['font.size'] = 16
# plt.bar(x,y,width=0.5,color='maroon')
# plt.ylabel('Count')
