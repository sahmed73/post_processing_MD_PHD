# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Mar  1 20:23:53 2025
"""

import magnolia.bondfile_parser as bfp


dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1"
filename = r"\bonds.reaxc"

atom_connectivity = bfp.parsebondfile(dirr+filename, mtypes=True)
#%%
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

neighbors = atom_connectivity['neighbours']
atypes    = atom_connectivity['atypes']
mtypes    = atom_connectivity['mtypes']
first_neighbors = neighbors[0]

def draw_backbone(molecule, carbons=[]):
    firststep_neigh     = list(neighbors.values())[0]
    firststep_graph     = nx.Graph(firststep_neigh)
    
    subgraph = firststep_graph.subgraph(molecule)
    H_nodes  = [x for x in molecule if atypes[x]==1]
    backbone = subgraph.copy()
    backbone.remove_nodes_from(H_nodes)
    
    ## Draw graph
    node_color = []
    for node in backbone.nodes:
        if node in carbons: node_color.append('black')
        else: node_color.append('gray')
        
    pos = nx.kamada_kawai_layout(backbone)
    plt.figure(figsize=[6,4])
    nx.draw(backbone,pos=pos,node_size=150,node_color=node_color,
            edgecolors='k')
    plt.show()

def find_first_ch_bond_break(atom_connectivity):
    neighbors = atom_connectivity['neighbours']  # Adjacency list per time step
    atypes = atom_connectivity['atypes']  # Atom ID to atom type mapping
    mtypes = atom_connectivity['mtypes']  # Atom ID to molecule ID mapping
    
    molecule_ch_bond_breaks = {}  # Store first break step for each molecule
    
    sorted_steps = sorted(neighbors.keys())  # Get sorted step keys
    
    for i in range(1, len(sorted_steps)):
        step = sorted_steps[i]
        prev_step = sorted_steps[i - 1]
        adjacency = neighbors[step]
        prev_adjacency = neighbors[prev_step]
        
        for atom_id, bonded_atoms in adjacency.items():
            if atypes[atom_id] == 2:  # Carbon (C) atom
                molecule_id = mtypes[atom_id]  # Get molecule ID
                
                # Count the number of Hydrogen (H) bonds
                initial_hydrogens = sum(1 for neighbor in bonded_atoms if atypes[neighbor] == 1)
                prev_bonded_atoms = prev_adjacency.get(atom_id, set())
                prev_hydrogens = sum(1 for neighbor in prev_bonded_atoms if atypes[neighbor] == 1)
                
                if prev_hydrogens > initial_hydrogens:
                    if molecule_id not in molecule_ch_bond_breaks:
                        molecule_ch_bond_breaks[molecule_id] = atom_id #(step*0.25/1000-325, atom_id)
    
    return molecule_ch_bond_breaks 
    

for step, neigh in neighbors.items():
    G = nx.Graph(neigh)
    molecules = bfp.get_molecules(neigh)
    
    for molecule in molecules:
        subgraph = G.subgraph(molecule)
        tert_C = []
        for atom in molecule:
            if atypes[atom] == 2:
                children = subgraph.neighbors(atom)
                C_children = [x for x in children if atypes[x]==2]
                N_Carbon = len(C_children)
                if N_Carbon==3:
                    tert_C.append(atom)
        if N_Carbon==3: draw_backbone(molecule, tert_C)
    break


carbon_types = {}  # Store classification of each carbon

for parent, children in first_neighbors.items():
    if atypes[parent] == 2:  # Only classify Carbon (C=2)
        bonded_carbons = [x for x in children if atypes[x] == 2] 
        num_carbons = len(bonded_carbons)

        if num_carbons == 1:
            carbon_types[parent] = "Primary"
        elif num_carbons == 2:
            carbon_types[parent] = "Secondary"
        elif num_carbons == 3:
            carbon_types[parent] = "Tertiary"

molecule_ch_bond_breaks = find_first_ch_bond_break(atom_connectivity)
#%%
degree = {}
result = {'Primary':0, 'Secondary':0, 'Tertiary':0}
for mid, aid in molecule_ch_bond_breaks.items():
    degree[aid]=carbon_types[aid]
    result[carbon_types[aid]]+=1
    
for ctype in result:
    result[ctype] = result[ctype]

fig, ax = plt.subplots(dpi=350)
ax.bar(result.keys(), result.values(), color='tab:red', edgecolor='k')
ax.set_ylabel('# first C-H bond breaking')