# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Mar 25 20:07:52 2024
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import networkx as nx
from networkx.algorithms import isomorphism
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

dirr      = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600"
sim_dirr  = ['Sim-1', 'Sim-2', 'Sim-3']

bonddata = {}
bo_cutoff = 0.3
for sim in sim_dirr:
    bondfilepath = dirr+'\\'+sim+'\\bonds.reaxc'
    bonddata[sim]=bfp.parsebondfile(bondfilepath, cutoff=bo_cutoff,
                                    bo=True)
#%%

## Collect specific molecule_graph
## Here I am only collecting 75 PAO from three sim, 25 from each
start = time.time()

mol_graphs  = []
ref = None

for sim in sim_dirr:
    bondfilepath = dirr+'\\'+sim+'\\bonds.reaxc'
    fsbd     = bfp.parsebondfile(bondfilepath,firststep=True) #only the first step
    fsneigh  = fsbd['neighbours'] # fsbd=first step bond data
    atypes   = fsbd['atypes']    
    fs_molecules= bfp.get_molecules(fsneigh)
    
    # get list of molecule_graph from the first steps
    molecule_length = 92
    fsneigh_graph = nx.Graph(fsneigh) 
    for molecule in fs_molecules:
        if len(molecule)==molecule_length:
            mol_graph = fsneigh_graph.subgraph(molecule)
            mol_graphs.append(mol_graph)
    break

ref = mol_graphs[0]            

def process_graph(mol_graph, ref, atypes):
    atom_map = {}
    site_map = {}
    matcher = isomorphism.GraphMatcher(mol_graph, ref)
    if matcher.is_isomorphic():
        mapping = matcher.mapping
        
        for node in mol_graph.nodes:
            atom_map[node] = mapping[node]
            
            # Assuming you have a way to get the neighbors
            # Corrected from 'neighbours' to 'neighbors' to match NetworkX's API
            if atypes[node] == 1:  # H
                connected_carbon = list(mol_graph.neighbors(node))[0]
                site_map[node] = mapping[connected_carbon]
            else:
                site_map[node] = mapping[node]
    return atom_map, site_map

# Assuming ref and mol_graphs are already defined

def parallel_process_graphs(mol_graphs, ref, atypes):
    combined_atom_map = {}
    combined_site_map = {}
    with ThreadPoolExecutor() as executor:
        futures = []
        for mol_graph in mol_graphs:
            futures.append(executor.submit(process_graph, mol_graph, ref, atypes))
        
        for future in as_completed(futures):
            atom_map, site_map = future.result()
            combined_atom_map.update(atom_map)
            combined_site_map.update(site_map)
    
    return combined_atom_map, combined_site_map

# Now call the function to process all graphs in parallel
atom_map, site_map = parallel_process_graphs(mol_graphs, ref, atypes)

end = time.time()
print(end-start)
#%%
