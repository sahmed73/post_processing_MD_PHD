# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Dec 19 00:51:03 2024
"""  

import numpy as np
from periodictable import elements
import networkx as nx
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import magnolia.MD_Converter as mdc
import os


# Example function to get covalent radius (can be replaced with a custom table)
def get_covalent_radius(symbol):
    return elements.symbol(symbol).covalent_radius

# Parse XYZ file
def read_xyz(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    num_atoms = int(lines[0].strip())
    data = [line.split() for line in lines[2:2 + num_atoms]]
    atoms = [line[0] for line in data]
    coordinates = np.array([[float(x) for x in line[1:]] for line in data])
    return atoms, coordinates

# Calculate bonds
def calculate_bonds(atoms, coordinates, tolerance=0.2):
    num_atoms = len(atoms)
    bonds = []
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            dist = np.linalg.norm(coordinates[i] - coordinates[j])
            r_i = get_covalent_radius(atoms[i])
            r_j = get_covalent_radius(atoms[j])
            if dist <= r_i + r_j + tolerance:
                bonds.append((i, j))  # Bond between atom i and j
    return bonds

# Main
dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\A0004'
xyzfile=dirr+r'\pre\A0004_optimized_pre.xyz'

which=1
os.makedirs(dirr+f'\\post{which}')
outfile=dirr+rf'\post{which}\A0004_un-optimized_post{which}.xyz'

symbols, coordinates = read_xyz(xyzfile)
bonds = calculate_bonds(symbols, coordinates)
    
##
G = nx.Graph()
G.add_edges_from(bonds)

phenolic_oxygens = []
rings = [ring for ring in nx.cycle_basis(G) if len(ring) == 6]
benzene_rings = [
    ring for ring in rings if all(symbols[node] == 'C' for node in ring)
]

# Step 3: Check for oxygens connected to benzene rings
for ring in benzene_rings:
    for node in ring:
        for neighbor in G.neighbors(node):
            if symbols[neighbor] == 'O':
                phenolic_oxygens.append(neighbor)

print(f'There are {len(phenolic_oxygens)} phenolic oxygens: {phenolic_oxygens}')                
                    
delete_atoms = []
for oxygen in phenolic_oxygens:
    for child in G.neighbors(oxygen):
        if symbols[child] == 'H':
            delete_atoms.append(child)
            break

updated_symbols = []
updated_coordinates = []
for i in range(len(coordinates)):
    if i != delete_atoms[which-1]:
        updated_symbols.append(symbols[i])
        updated_coordinates.append(coordinates[i])
mdc.write_xyzfile(outfile, updated_symbols, updated_coordinates)

#%% dreaw molecules
fig , ax=plt.subplots(1,2,dpi=350)

# pre
types = [{'H':1, 'C':2, 'O':3}[x] for x in symbols]
atom_types = dict(zip(range(0,len(coordinates)), types))
man.draw_rdkit2D_from_graph(G, atom_types, "HCO",ax=ax[0],
                            atom_label=False,
                            highlight_atoms=[phenolic_oxygens[which-1]])

ax[0].set_title('Pre',fontsize=12)

# post
updated_G = G.copy()
updated_G.remove_node(delete_atoms[which-1])
types = [{'H':1, 'C':2, 'O':3}[x] for x in updated_symbols]
atom_types = dict(zip(range(0,len(updated_coordinates)), types))
man.draw_rdkit2D_from_graph(updated_G, atom_types, "HCO",ax=ax[1],
                            atom_label=False,
                            highlight_atoms=[phenolic_oxygens[which-1]])
ax[1].set_title(f'Post_{which}',fontsize=12)
