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
import seaborn as sns

baseoil = "PAO4"
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(baseoil,baseoil)
bondfilepath = directory+'\\bonds.reaxc'
cutoff = 0.35
bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%% seeking from the last step
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
atomsymbols = 'HCO'

seek = 'H6C3O2'
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
#%% carbon skeleton
def draw_backbone(molecule_type,frame=0,carbons=[],title=''):
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
#%% Carbon Sources
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\pathway'
for i, track in enumerate(seek_molecules):
    carbons = [x for x in track if atypes[x]==2]
    molecule_type = mtypes[carbons[0]]
    draw_backbone(molecule_type,carbons=carbons)
    plt.title(f'{i}:  {seek}',fontsize=30)
    plt.savefig(savedir+f'\\backbone_{i}.png', dpi=500, bbox_inches='tight')