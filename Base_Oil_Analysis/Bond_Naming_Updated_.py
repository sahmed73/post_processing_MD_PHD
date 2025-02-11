# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Dec 19 22:38:15 2023
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
from networkx.algorithms import isomorphism
from networkx.drawing.nx_agraph import graphviz_layout
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import pandas as pd
import seaborn as sns
import numpy as np

def map_isomorphic_nodes_to_unique_ids(graphs,reference_graph):
    """
    Maps each node in a list of isomorphic graphs to a unique ID based on node position.

    :param graphs: List of NetworkX graphs that are isomorphic to each other.
    :return: Dictionary with node IDs as keys and unique IDs (1-30) as values.
    """
    if not graphs:
        return {}

    node_to_unique_id = {node: unique_id for unique_id, node in enumerate(reference_graph.nodes, start=1)}

    # Map nodes of the other graphs
    for graph in graphs[1:]:
        iso_matcher = isomorphism.GraphMatcher(reference_graph, graph)
        for iso_map in iso_matcher.isomorphisms_iter():
            # Apply the unique ID from the reference node to the corresponding node in this graph
            for ref_node, node in iso_map.items():
                node_to_unique_id[node] = node_to_unique_id[ref_node]
            break  # Consider only the first mapping

    return node_to_unique_id
#%%
sim_dir = ['Sim-1','Sim-2','Sim-3']
baseoil = "PAO4"
ref_flag = True
count = {}

for sim in sim_dir:
    dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(baseoil,baseoil,sim)
    
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\{}'.format(sim)
    
    bondfilepath = dirr+'\\bonds.reaxc'
    cutoff = 0.45
    bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, bo=True)
    
    neighbors = bonddata['neighbours']
    atypes    = bonddata['atypes']
    bondorders= bonddata['bondorders']
    atomsymbols = 'HCO'
    
    
    ## getting a list of baseoil
    oils = []
    firrstneigh = neighbors[0]
    main_graph = nx.Graph(firrstneigh)
    molecules = bfp.get_molecules(firrstneigh)
    for molecule in molecules:
        molGraph = main_graph.subgraph(molecule).copy()
        if len(molecule)==92:
            nodes_to_remove=[node for node in molGraph.nodes if atypes[node]==1]
            molGraph.remove_nodes_from(nodes_to_remove)
            oils.append(molGraph)
    
    if ref_flag:
        reference_graph = oils[0].copy()
        ref_flag = False
    
    carbon_enum = map_isomorphic_nodes_to_unique_ids(oils,reference_graph)
    for ii, values in enumerate(oils):
        oil = values.copy()
        labels = {node:carbon_enum[node] for node in oil.nodes}
        node_color = ['red' if labels[node]==11 else 'lightblue' for node in oil]
        pos = graphviz_layout(oil, prog='neato')
        nx.draw(oil,with_labels=True,labels=labels,node_color=node_color,pos=pos)
        plt.title(f'Carbon Skeleton: {ii+1}')
        plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs\carbon_skeleton_{}_{}'.format(ii+1,sim),dpi=500,bbox_inches='tight')
        plt.show()
    
    if not count:
        for u,v in oil.edges:
            uu = carbon_enum[u]
            vv = carbon_enum[v]
            pair = tuple(sorted((uu,vv)))
            count[pair]=0
    
    bo_cutoff = 0.45
    for step in bondorders:
        for atom1 in bondorders[step]:
            if atypes[atom1]!=2: continue
            for atom2 in bondorders[step][atom1]:
                if atypes[atom2]!=2: continue
                bond_order = bondorders[step][atom1].get(atom2,0)
                if bond_order<bo_cutoff:
                    uu = carbon_enum[atom1]
                    vv = carbon_enum[atom2]
                    pair = tuple(sorted((uu,vv)))
                    if pair in count:
                        count[pair]+=1
#%%
special= 4

sorted_items = sorted(count.keys(), key=lambda x: count[x], reverse=True)
top_counts = list(sorted_items)[:special]
colors = ['red' if k in top_counts else 'tab:blue' for k in count]
    
plt.bar(list(map(str,count.keys())),count.values(),color=colors)
plt.xticks(rotation=90,fontsize=7)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\carbon_bond_break_count',dpi=500,bbox_inches='tight')
plt.show()

top_nodes = list({e for tup in top_counts for e in tup})

labels = {node:carbon_enum[node] for node in oil.nodes}
node_color = ['red' if labels[node] in top_nodes else 'tab:blue' for node in oil]
edge_colors = ['red' if tuple(sorted((labels[x],labels[y]))) in top_counts else 'k' for x,y in oil.edges]
edge_widths = [3 if tuple(sorted((labels[x],labels[y]))) in top_counts else 1 for x,y in oil.edges]

plt.figure(figsize=[10,6])
pos = graphviz_layout(oil, prog='neato')
nx.draw(oil,with_labels=True,labels=labels,node_color=node_color,pos=pos)
nx.draw_networkx_edges(oil, pos, edge_color=edge_colors,width=edge_widths)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs\specieal_carbon_skeleton_{}'.format(baseoil),dpi=500,bbox_inches='tight')
plt.show()