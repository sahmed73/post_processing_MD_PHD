# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Jan 29 02:10:58 2024
"""
import networkx as nx
import magnolia.bondfile_parser as bfp

def draw_backbone(bonddata,molecule,frame=0,carbons=[],title=''):
    neighbors = bonddata['neighbours']
    atypes    = bonddata['atypes']
    mtypes    = bonddata['mtypes']
    bondorders= bonddata['bondorders']
    
    if isinstance(molecule_type, int):
        firststep_neigh     = list(neighbors.values())[0]
        firststep_graph     = nx.Graph(firststep_neigh)
        molecule            = [k for k,v in mtypes.items() if v==molecule_type]
    else:
        firststep_neigh     = list(neighbors.values())[frame]
        firststep_graph     = nx.Graph(firststep_neigh)
        molecule            = molecule_type.copy()
    
    subgraph = firststep_graph.subgraph(molecule)
    H_nodes  = [x for x in molecule if atypes[x]==1]
    backbone = subgraph.copy()
    backbone.remove_nodes_from(H_nodes)
    center   = nx.center(backbone)
    print(f'center = {center}')
    
    ## Draw graph
    node_color = []
    for node in backbone.nodes:
        if node in carbons: node_color.append('gold')
        elif atypes[node]==3: node_color.append('red')
        else: node_color.append('tab:blue')
        
    pos = nx.kamada_kawai_layout(backbone)
    plt.figure(figsize=[6,4])
    nx.draw(backbone,pos=pos,node_size=150,node_color=node_color)
    plt.title(title,fontsize=30)