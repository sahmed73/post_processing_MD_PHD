# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Jan 25 12:20:35 2024
"""
## update species pathway tracker

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import numpy as np
import random
import os

baseoil = "PAO4"
cutoff = 0.30
bonddata = {}
simdir = ['Sim-1']#,'Sim-2','Sim-3']
for sim in simdir:
    directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(baseoil,baseoil,sim)
    bondfilepath = directory+'\\bonds.reaxc'
    bonddata[sim] = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%% seeking from the last step
sim = 'Sim-1'

neighbors = bonddata[sim]['neighbours']
atypes    = bonddata[sim]['atypes']
mtypes    = bonddata[sim]['mtypes']
bondorders= bonddata[sim]['bondorders']
atomsymbols = 'HCO'

seek = 'H2CO'
seek_molecules = set()

laststep_neigh     = list(neighbors.values())[-1]
laststep           = list(neighbors.keys())[-1]
laststep_molecules = list(bfp.get_molecules(laststep_neigh))

for molecule in laststep_molecules:
    species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
    if species==seek and frozenset(molecule) not in seek_molecules:
        seek_molecules.add(frozenset(molecule))
        print(molecule)
        for atom1 in molecule:
            for atom2 in molecule:
                if atom1>atom2:
                    sym1 = atomsymbols[atypes[atom1]-1]
                    sym2 = atomsymbols[atypes[atom2]-1]
                    bb = bondorders[laststep][atom1].get(atom2,None)
                    if bb: print((sym1,sym2),bb)
#%% tracking a single seeking molecule in reveres step order
def find_all_largest_consecutive_sequences(arr):
    if not arr:
        return []

    arr.sort()  # Sort the input list

    longest_sequences = []  # Initialize a list to store all longest sequences
    current_sequence = [arr[0]]  # Initialize the current sequence with the first element

    for i in range(1, len(arr)):
        if arr[i] == arr[i - 1] + 1:
            current_sequence.append(arr[i])
        else:
            # Check if the current sequence is longer than or equal to the longest sequence(s)
            if not longest_sequences or len(current_sequence) == len(longest_sequences[0]):
                longest_sequences.append(current_sequence.copy())
            elif len(current_sequence) > len(longest_sequences[0]):
                longest_sequences = [current_sequence.copy()]
            current_sequence = [arr[i]]

    # Check again after the loop ends in case the longest sequence(s) is/are at the end
    if not longest_sequences or len(current_sequence) == len(longest_sequences[0]):
        longest_sequences.append(current_sequence.copy())
    elif len(current_sequence) > len(longest_sequences[0]):
        longest_sequences = [current_sequence.copy()]

    return np.array(longest_sequences)

## only tracking 1
print(seek)
for i in range(len(seek_molecules)):
    if i!=10: continue
    track          = list(seek_molecules)[i]
    exist_in_frame = []
    
    reversed_keys = list(neighbors.keys())[::-1]
    total_frame   = len(reversed_keys)
    for f, key in enumerate(reversed_keys):
        frame = total_frame-f-1
        neigh = neighbors[key]
        molecules = bfp.get_molecules(neigh)
        for molecule in molecules:
            if track == molecule:
                exist_in_frame.append(frame)
                bingo = True
            
                
    result = find_all_largest_consecutive_sequences(exist_in_frame)
    appeared, life = result.shape
    print(f"seeking id: {i}: appeared = {appeared}, life = {life}")
#%%
# working with single track
parents = {}
for frame, (step,neigh) in enumerate(neighbors.items()):
    parents[frame]=[]
    molecules = bfp.get_molecules(neigh)
    for molecule in molecules:
        if track & molecule:
            parents[frame].append(molecule)
#%%
#draw molecule
from rdkit import Chem
db_cutoff  = 1.40 # double bond cutoff

def nx_to_rdkit(G,step):
    mol = Chem.RWMol()  # Create an editable RDKit molecule
    
    # Add atoms to the molecule
    atom_indices = {}
    atomic_number = [1,6,8]
    for node in G.nodes:
        atom_num = atomic_number[atypes[node]-1]
        atom = Chem.Atom(atom_num)
        idx = mol.AddAtom(atom)
        atom_indices[node] = idx
    
    # Add bonds to the molecule
    for edge in G.edges():
        start, end = edge
        bo = bondorders[step][start].get(end,0)
        if bo>db_cutoff:
            mol.AddBond(atom_indices[start], atom_indices[end], Chem.BondType.DOUBLE)
        else:
            mol.AddBond(atom_indices[start], atom_indices[end], Chem.BondType.SINGLE)
    
    # Convert the editable molecule to a finalized molecule
    mol = mol.GetMol()
    
    return mol

def draw_backbone(molecule,frame,show_H=False):
    species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
    step = list(neighbors.keys())[frame]
    neigh = neighbors[step]
    G     = nx.Graph(neigh)    
    subgraph = G.subgraph(molecule)
    H_nodes  = [x for x in molecule if atypes[x]==1]
    backbone = subgraph.copy()
        
    if not show_H:
        backbone.remove_nodes_from(H_nodes)
    
    ## smiles string
    rdkit_mol = nx_to_rdkit(subgraph,step)
    smiles = Chem.MolToSmiles(rdkit_mol)
        
    ## Draw graph
    node_color = []
    node_size  = []
    node_shape = []
    for node in backbone.nodes:
        if node in track:
            # node_color.append('gold')
            if atypes[node]==3:
                node_size.append(100)
                node_color.append('tab:red')
                node_shape.append('blue')
            elif atypes[node]==2:
                node_color.append('tab:grey')
                node_size.append(180)
                node_shape.append('blue')
            else:
                node_size.append(30)
                node_color.append('white')
                node_shape.append('blue')
        elif atypes[node]==3:
            node_color.append('tab:red')
            node_size.append(100)
            node_shape.append('k')
        elif atypes[node]==2:
            node_color.append('tab:grey')
            node_size.append(180)
            node_shape.append('k')
        elif atypes[node]==1:
            node_color.append('white')
            node_size.append(30)
            node_shape.append('k')
        else: 
            node_color.append('gold')
            node_size.append(350)
            node_shape.append('k')
    
    edge_color = []
    edge_thickness = []
    for u,v in backbone.edges:
        if bondorders[step][u][v]>db_cutoff:
            edge_color.append('tab:red')
            edge_thickness.append(3)
        else:
            edge_color.append('k')
            edge_thickness.append(1.0)
    
    pos = nx.kamada_kawai_layout(backbone)
    plt.figure(figsize=[10,8])
    
    nx.draw(backbone,pos=pos,node_size=node_size,node_color=node_color,
            edgecolors=node_shape, edge_color=edge_color,
            width=edge_thickness)
    plt.title(f'frame-{frame}: {species}',fontsize=20)
    return smiles    
#%% Analyze the parents
save = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs'
spec = bfp.get_molecular_formula(track, atypes, atomsymbols)
dirr = save+f'\\4_20241303_{baseoil}_{spec}'
if not os.path.exists(dirr):
    os.makedirs(dirr)

check_list = []    
for frame in range(5301):
    if frame%100==0: print(frame)
    parent = parents[frame]
    for molecule in parent:
        species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
        if species not in check_list:
            smiles = draw_backbone(molecule, frame,show_H=True)
            print(frame,smiles)
            plt.savefig(dirr+f'\\graph_{frame}_{species}.png',dpi=500,bbox_inches='tight')
            plt.close()
            check_list.append(species)

#%% Only Carbon Common
save = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs'
spec = bfp.get_molecular_formula(track, atypes, atomsymbols)
dirr = save+f'\\{baseoil}_{spec}_carbon_based'
if not os.path.exists(dirr):
    os.makedirs(dirr)
    
check_list = []
for frame in range(5301):
    if frame%100==0: print(frame)
    parent = parents[frame]
    for molecule in parent:
        common = track & molecule
        types  = [atypes[x] for x in common]
        species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
        if species not in check_list and 2 in types:
            smiles = draw_backbone(molecule, frame,show_H=True)
            print(frame,smiles)
            plt.savefig(dirr+f'\\graph_{frame}_{species}.png',dpi=300,bbox_inches='tight')
            plt.close()
            check_list.append(species)    