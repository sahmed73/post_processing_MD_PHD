# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Apr  5 09:43:54 2024
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import networkx as nx
from networkx.algorithms import isomorphism
from concurrent.futures import ProcessPoolExecutor

dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600"

nsim = 6
sim_dirr = [f'Sim-{x+1}' for x in range(nsim)]

# Define a function for concurrent execution
def parse_bond_file(sim):
    bondfilepath = dirr + '\\' + sim + '\\bonds.reaxc'
    return bfp.parsebondfile(bondfilepath, cutoff=0.3, bo=True, mtypes=True)

# Use ProcessPoolExecutor to parallelize bond file parsing
with ProcessPoolExecutor() as executor:
    results = list(executor.map(parse_bond_file, sim_dirr))

# Convert results into a dictionary
bonddata = dict(zip(sim_dirr, results))

#%%
## creating atom_map and site_map

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
            

atom_map = {}
site_map = {}
if ref is None:
    ref=mol_graphs[0]

for i, mol_graph in enumerate(mol_graphs):
    print(i)
    matcher = isomorphism.GraphMatcher(mol_graph, ref)
    if matcher.is_isomorphic():
        mapping = matcher.mapping
        
        for node in mol_graph.nodes:
            atom_map[node] = mapping[node]
            
            # site mapping
            if atypes[node] == 1: # H
                connected_carbon = list(mol_graph.neighbors(node))[0]
                site_map[node] = mapping[connected_carbon]
            else:
                site_map[node] = mapping[node]

end = time.time()
print(end-start)
#%%
start = time.time()
site_time = {}
count = 1
data = []
m_sim_type_check = []
for sim in sim_dirr:
    print(sim)
    neighbors = bonddata[sim]['neighbours']
    atypes    = bonddata[sim]['atypes']
    mtypes    = bonddata[sim]['mtypes']
    
    seek = 'H61C30'
    for step, neigh in neighbors.items():
        molecules = bfp.get_molecules(neigh)
        for molecule in molecules:
            species = bfp.get_molecular_formula(molecule, atypes, 'HCO')
            m_sim_type = sim+'-'+str(mtypes[list(molecule)[0]])
            if species==seek and m_sim_type_check.count(m_sim_type)<=count:
                data.append((molecule,step))
                m_sim_type_check.append(m_sim_type)
end = time.time()
print(end-start)
#%%
first = []
for molecule,step in data:
    fsbd     = bfp.parsebondfile(bondfilepath,firststep=True) #only the first step
    fsneigh  = fsbd['neighbours'] # fsbd=first step bond data
    atypes   = fsbd['atypes']    
    fs_molecules= bfp.get_molecules(fsneigh)
    
    # get list of molecule_graph from the first steps
    molecule_length = 92
    fsneigh_graph = nx.Graph(fsneigh)
    PAOs = []
    for fs_molecule in fs_molecules:
        if len(fs_molecule)==molecule_length:
            PAOs.append(fs_molecule)
    
    for PAO in PAOs:
        uncommon = PAO^molecule
        if len(uncommon)==1:
            single = list(uncommon)[0]
            first.append(site_map[single])
 
plt.rc('font', size=14) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=18) # Axes label size
plt.rc('xtick', labelsize=16) # X-axis tick label size
plt.rc('ytick', labelsize=16) # Y-axis tick label size
plt.rc('legend', fontsize=12) # Legend fontsize

first = np.array(first)
unique_elements, counts = np.unique(first, return_counts=True)

fig, ax = plt.subplots()
bars = ax.bar(unique_elements-min(unique_elements)+1,counts,
        color='tab:red',edgecolor='black')
ax.set_xlabel('H site')
ax.set_ylabel('Count')
ax.set_xticks(range(0,31,5))
# ax.set_ylim(top=6)
# ax.set_yticks(range(6))

# Adding a star sign on top of the highest bars
# max_height = max(counts)
# for bar in bars:
#     height = bar.get_height()
#     if height == max_height:
#         ax.text(bar.get_x() + bar.get_width() / 2., height, 
#                 '$\star$', ha='center', va='bottom', fontsize=30, color='red')
        
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\H_breaking', dpi=300,
            bbox_inches='tight')