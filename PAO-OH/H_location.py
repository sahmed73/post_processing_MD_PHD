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
import magnolia.access_ucm_cluster as ucm

AO = 1
directory = r'/mnt/borgstore/amartini/sahmed73/ARLOI-V1/4_ARLOI_24-09-07--22-34-03/20_PAO-OH_15_A0001/Production/Sim-1'
filename  = '/bonds.reaxc'
bondfile  = ucm.local_copy(directory+filename)
atomsymbols = ['H','C','O']
atominfo = bfp.parsebondfile(bondfile,bo=True)
#%%
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
                    
#%%            
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

y = [len(positive), len(negative)]
x = ['Hydroxyl-H', 'Nonhydroxyl-H']

plt.subplots(dpi=350)
plt.rcParams['font.size'] = 16
plt.bar(x,y,width=0.5,color='maroon')
plt.ylabel('Count')
plt.title(f'A000{AO}\n')