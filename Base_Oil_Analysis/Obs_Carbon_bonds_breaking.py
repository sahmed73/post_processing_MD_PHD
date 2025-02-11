# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Dec  3 13:23:57 2023
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import pandas as pd
import seaborn as sns
import numpy as np

baseoil = "PAO4"
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(baseoil,baseoil)
bondfilepath = directory+'\\bonds.reaxc'
cutoff = 0.45
bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
atomsymbols = 'HCO'

## bfs helper function
def assign_ids_by_traversal(tree, source, ref):
    unique_id = 1
    node_id_mapping = {}

    for node in nx.bfs_tree(tree, source=source).nodes():
        if node not in node_id_mapping:
            node_id_mapping[node] = unique_id
            unique_id += 1
    return node_id_mapping


## getting a list of baseoil
oils = []
firrstneigh = neighbors[0]
main_graph = nx.Graph(firrstneigh)
molecules = bfp.get_molecules(firrstneigh)
for molecule in molecules:
    molGraph = main_graph.subgraph(molecule)
    if len(molecule)==92:
        oils.append(molGraph)

## number all the carbon using bfs
carbon_enum = {}
first_flag = True

for ii, values in enumerate(oils):
    oil = values.copy()
    nodes_to_remove=[node for node in oil.nodes if atypes[node]==1]
    oil.remove_nodes_from(nodes_to_remove)
    
    if first_flag:
        first_flag=False
        ref = oil.copy()
    
    center = nx.center(oil)[0] # take single
    carbon_enum |= assign_ids_by_traversal(oil, center)
    labels = {node:carbon_enum[node] for node in oil.nodes}
    node_color = ['red' if node==center else 'lightgrey' for node in oil.nodes]
    nx.draw_spring(oil,with_labels=True,node_color=node_color,labels=labels)
    plt.title(f'Carbon Skeleton: {ii+1}')
    plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs\carbon_skeleton_{}'.format(ii+1),dpi=500,bbox_inches='tight')
    plt.show()
#%%
bonds = []
for u,v in ref.edges:
    uu = carbon_enum[u]
    vv = carbon_enum[v]
    bond = tuple(sorted((uu,vv)))
    bonds.append(bond)
    
counts = {}
for step in bondorders:
    for atom1 in bondorders[step]:
        if atypes[atom1]!=2: continue
        for atom2 in bondorders[step][atom1]:
            if atypes[atom2]!=2: continue
            bond_order = bondorders[step][atom1].get(atom2,0)
            cnum1 = carbon_enum[atom1]
            cnum2 = carbon_enum[atom2]
            if bond_order<cutoff:
                pair = tuple(sorted((cnum1,cnum2)))
                if pair not in bonds: continue
                if pair in counts:
                    counts[pair]+=1
                else:
                    counts[pair]=1
                    
#%%
center = nx.center(ref)[0] # take single
labels = {node:carbon_enum[node] for node in ref.nodes}
node_color = ['red' if node==center else 'lightgrey' for node in ref.nodes]
pos = nx.spring_layout(ref)
# nx.draw(ref,with_labels=True,node_color=node_color,labels=labels,pos=pos)
# edge_labels = nx.get_edge_attributes(ref, 'count')
# nx.draw_networkx_edge_labels(ref, pos, edge_labels=edge_labels)
# plt.show()
aa = list((map(str,counts.keys())))
plt.bar(aa,counts.values())
# plt.xticks(aa,rotation=90)
# plt.ylim(0,10)
# plt.show()
#%% Draw a single one
plt.figure(figsize=[10,10])
labels = {node:carbon_enum[node] for node in oil.nodes}
node_color = ['red' if node==center else 'lightgrey' for node in oil.nodes]
nx.draw_spring(oil,with_labels=True,node_color=node_color,labels=labels)
plt.show()