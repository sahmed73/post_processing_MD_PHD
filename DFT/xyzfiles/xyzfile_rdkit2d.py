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
dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0003\DataFile\Single_Molecules'
xyzfile=dirr+r'\A0003_pre.xyz'

atoms, coordinates = read_xyz(xyzfile)
bonds = calculate_bonds(atoms, coordinates)

# Output bonds
print("Bonds (atom indices):")
for bond in bonds:
    print(bond)
    
    
##
G = nx.Graph()
G.add_edges_from(bonds)
fig,ax=plt.subplots(dpi=1200)
types = [{'H':1, 'C':2, 'O':3}[x] for x in atoms]
atom_types = dict(zip(range(0,len(coordinates)), types))
man.draw_rdkit2D_from_graph(G, atom_types, "HCO",ax=ax,
                            atom_label=False, highlight_atoms=None)
fig.savefig("2d_molecule.png",dpi=1200)
        