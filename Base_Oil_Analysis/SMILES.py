# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov  6 21:55:13 2023
"""

import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
import random
import numpy as np
import networkx as nx
from rdkit import Chem

### Directory ###
base_oil = 'PAO4'
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(base_oil,base_oil)
filename  = "\\bonds.reaxc"
bondfilepath = directory+filename

### Parsing Bondfile ###
bonddata = bfp.parsebondfile(bondfilepath,bo=True)
#%%
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
# bondorders = bonddata['bondorders']
atomsymbols= ['H','C','O']
atomic_num   = [1,6,8]

def drawGraph(G):
    # Layout for the nodes in the graph
    pos = nx.spring_layout(G)  

    # Drawing the graph
    nx.draw_networkx_nodes(G, pos, node_size=700)
    nx.draw_networkx_edges(G, pos, width=6)
    nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")

    # Draw edge labels (weights)
    edge_labels = nx.get_edge_attributes(G, "weight")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    # Turn off axis and show the graph
    plt.axis("off")
    plt.show()

images = []

smiles_array = []
for step, neigh in neighbours.items():
    graph     = nx.Graph(neigh)
    molecules = nx.connected_components(graph)
    bo        = bondorders[step]
    bo_array, labels = bfp.get_bondtypes(bo, n_clusters=2)
    for molecule in molecules:
        species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
        if species!='H2CO2':
            continue
        rdkitMol = Chem.RWMol()
        subgraph = graph.subgraph(molecule)
        for u,v in subgraph.edges():
            bo_order = bo[u][v]
            if bo_order in bo_array[labels==1]:
                bond_type = 1
            elif bo_order in bo_array[labels==2]:
                bond_type = 2
            elif bo_order in bo_array[labels==3]:
                bond_type = 3
            else:
                print('bo_order out of range!!')
            subgraph[u][v]['bond_type']=bond_type
        smiles   = bfp.atomConnectivity2smiles(subgraph,atypes,atomic_num)
        
        
        molecule = Chem.MolFromSmiles(smiles)
        # Draw the molecule
        img = Chem.Draw.MolToImage(molecule)
        Chem.Draw.MolToFile(molecule, 'molecule.png')
        images.append(img)
        smiles_array.append(smiles)
        print(smiles)
        