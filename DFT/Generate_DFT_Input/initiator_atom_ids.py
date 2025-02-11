# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Dec 25 01:10:39 2024
"""

## for only PAO radical and antioxidant system

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man

# Parse LAMMPS data file
def read_lammps_data(file_path):
    symbols = []
    coordinates = []
    bonds = []
    molecule_ids = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        atoms_section = False
        bonds_section = False
        for line in lines:
            if "Atoms" in line:
                atoms_section = True
                continue
            if "Bonds" in line:
                atoms_section = False
                bonds_section = True
                continue
            if atoms_section and len(line.split()) > 0:
                # Example format: atom_id molecule_id atom_type x y z
                data = line.split()
                x, y, z = map(float, data[3:6])
                molecule_ids.append(int(data[1]))
                symbol = line.strip().split('/')[-1]
                symbols.append(symbol)
                coordinates.append([x, y, z])
            if bonds_section and len(line.split()) > 0:
                # Example format: bond_id bond_type atom1 atom2\
                if "Angles" in line:
                    bonds_section=False
                    break
                data = line.split()
                atom1, atom2 = map(int, data[2:4])
                bonds.append((atom1, atom2))
    return symbols, molecule_ids, np.array(coordinates), bonds

# Main logic
dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0004\DataFile\from_LUNAR\from_all2lmp'
pre_reaction_file = dirr + r'\pre_reaction_typed_IFF.data'

# 1 based indexing .. please be aware

symbols, molecule_ids, coordinates, bonds = read_lammps_data(pre_reaction_file)

# Create graph from bonds
G = nx.Graph()
G.add_edges_from(bonds)

# Identify phenolic oxygens
phenolic_oxygens = []
rings = [ring for ring in nx.cycle_basis(G) if len(ring) == 6]
benzene_rings = [
    ring for ring in rings if all(symbols[node-1] == 'C' for node in ring)
]

# Check for oxygens connected to benzene rings
for ring in benzene_rings:
    for node in ring:
        for neighbor in G.neighbors(node):
            if symbols[neighbor-1] == 'O':
                phenolic_oxygens.append(neighbor)

                
# add hydrogen attached to phenolic oxygen
initiator_atoms = []
for oxygen in phenolic_oxygens:
    for child in G.neighbors(oxygen):
        if symbols[child-1] == 'H':
            print(f'Phenolic hydroxyl hydrogen: {child}') 
            initiator_atoms.append(child)
            break

# identify the PAO radical Oxygen atom
for node in G:
    if molecule_ids[node-1]==1 and symbols[node-1]=='O':
        print(f'PAO radical oxygen: {node}')
        initiator_atoms.append(node)

print(f'initiator atoms: {initiator_atoms}')
#%% Draw molecules
fig, ax = plt.subplots(dpi=350)

# Pre
types = [{'H': 1, 'C': 2, 'O': 3}[x] for x in symbols]
atom_types = dict(zip(range(1, 1+len(coordinates)), types))
man.draw_rdkit2D_from_graph(G, atom_types, "HCO", ax=ax,
                            atom_label=False,
                            highlight_atoms=initiator_atoms)
#%% change in 

out_dir=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0004\DataFile\from_LUNAR\from_bond_react_merge'
rxnfile = out_dir+'\\pre1-post1_rxn-map_uncommented.txt'

# Read the file and modify the content
with open(rxnfile, 'r') as file:
    lines = file.readlines()

# Replace 'ID1' and 'ID2' with values from initiator_atoms
updated_lines = []
for line in lines:
    if 'ID1' in line or 'ID2' in line:
        line = line.replace('ID1', str(initiator_atoms[0]))
        line = line.replace('ID2', str(initiator_atoms[1]))
    updated_lines.append(line)

# Write the updated content back to the file
with open(rxnfile, 'w') as file:
    file.writelines(updated_lines)

print(f"Replaced 'ID1' and 'ID2' with {initiator_atoms} in the file.")
