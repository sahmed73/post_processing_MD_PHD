# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Feb  5 13:54:55 2025
"""

import numpy as np
from scipy.spatial import distance_matrix
from periodictable import elements
import networkx as nx
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man

# Function to get covalent radius
def get_covalent_radius(symbol):
    return elements.symbol(symbol).covalent_radius


# Function to calculate bonds
def calculate_bonds(atomic_info, tolerance=0.2):
    """
    Calculates bonds based on atomic coordinates and covalent radii.
    
    :param atomic_info: Dictionary {atom_id: [symbol, coordinates]} (0-based indices).
    :param tolerance: Extra distance allowance for bond formation.
    :return: List of bonded atom pairs (0-based indices).
    """
    atom_ids = list(atomic_info.keys())  # Get atom IDs
    num_atoms = len(atom_ids)
    bonds = []
    
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            atom_id_i = atom_ids[i]
            atom_id_j = atom_ids[j]

            symbol_i, coord_i = atomic_info[atom_id_i]
            symbol_j, coord_j = atomic_info[atom_id_j]

            dist = np.linalg.norm(np.array(coord_i) - np.array(coord_j))  # Euclidean distance
            r_i = get_covalent_radius(symbol_i)
            r_j = get_covalent_radius(symbol_j)

            if dist <= (r_i + r_j + tolerance):
                bonds.append((atom_id_i, atom_id_j))  # Store bond as 0-based indices

    return bonds


# Function to read XYZ file and store atomic info with 0-based indices
def read_xyz(filename):
    """Reads an XYZ file and returns a dictionary with 0-based atom IDs, symbols, and coordinates."""
    
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    num_atoms = int(lines[0].strip())
    atomic_info = {}  # Dictionary to store atom ID -> [Symbol, Coordinates]
    
    for i, line in enumerate(lines[2:num_atoms + 2]):  # Start from 0-based index
        parts = line.split()
        atomic_symbol = parts[0]
        coordinates = [float(parts[1]), float(parts[2]), float(parts[3])]
        
        atomic_info[i] = [atomic_symbol, coordinates]  # Store in dictionary

    return atomic_info


# Function to calculate gamma (steric hindrance factor)
def calculate_gamma(pao_file, ao_file, pao_atom_idx, ao_atom_idx, r_bond=1):
    """
    Computes the steric hindrance factor (gamma), where min_r is the smallest sphere that includes 
    all tert-butyl H, and n counts all H inside that sphere but excludes those directly bonded to the radical site.

    :param pao_file: Path to the PAO file (XYZ format).
    :param ao_file: Path to the AO file (XYZ format).
    :param pao_atom_idx: 0-based index of the reactive PAO atom.
    :param ao_atom_idx: 0-based index of the reactive AO atom.
    :param r_bond: Bond radius factor for steric hindrance.
    :return: Tuple (gamma, n1_neighbors, n2_neighbors).
    """
    atomic_info_pao = read_xyz(pao_file)
    atomic_info_ao = read_xyz(ao_file)
    
    bonds_pao = calculate_bonds(atomic_info_pao)
    bonds_ao = calculate_bonds(atomic_info_ao)
    
    # Compute pairwise distances within each molecule
    coords_pao = np.array([info[1] for info in atomic_info_pao.values()])
    coords_ao = np.array([info[1] for info in atomic_info_ao.values()])
    dist_pao = distance_matrix(coords_pao, coords_pao)
    dist_ao  = distance_matrix(coords_ao, coords_ao)
    
    # Find directly bonded atoms to the reaction site (exclude these H later)
    directly_bonded_pao = {j for i, j in bonds_pao if i == pao_atom_idx or j == pao_atom_idx}
    directly_bonded_ao = {j for i, j in bonds_ao if i == ao_atom_idx or j == ao_atom_idx}
    
    # Find tert-butyl hydrogen atoms
    tert_butyl_h_pao = find_tert_butyl_hydrogens(atomic_info_pao, bonds_pao)
    tert_butyl_h_ao = find_tert_butyl_hydrogens(atomic_info_ao, bonds_ao)
    
    # Step 1: Calculate min_r that includes all tert-butyl H
    if tert_butyl_h_pao:
        min_r_pao = max(dist_pao[pao_atom_idx, h] for h in tert_butyl_h_pao)
    else:
        min_r_pao = 0  # If no tert-butyl H exists, set to 0
    
    if tert_butyl_h_ao:
        min_r_ao = max(dist_ao[ao_atom_idx, h] for h in tert_butyl_h_ao)
    else:
        min_r_ao = 0  # If no tert-butyl H exists, set to 0
    
    # Step 2: Count all hydrogens inside min_r but exclude those directly bonded to the radical site
    n1_neighbors = [
        i for i, (symbol, _) in atomic_info_pao.items()
        if symbol == 'H' and dist_pao[pao_atom_idx, i] <= min_r_pao and i not in directly_bonded_pao
    ]
    
    n2_neighbors = [
        i for i, (symbol, _) in atomic_info_ao.items()
        if symbol == 'H' and dist_ao[ao_atom_idx, i] <= min_r_ao and i not in directly_bonded_ao
    ]

    # Count the number of hydrogen atoms within min_r
    n1 = len(n1_neighbors)
    n2 = len(n2_neighbors)
    
    # Compute steric hindrance factor
    print(f"r_bond: {r_bond}")
    print(f"n1: {n1}")
    print(f"n2: {n2}")
    print(f"min_r_pao: {min_r_pao}")
    print(f"min_r_ao: {min_r_ao}")
    gamma = (r_bond * (n1 + n2)) / (min_r_pao + min_r_ao) if (min_r_pao + min_r_ao) > 0 else 0
    
    return gamma, n1_neighbors, n2_neighbors


# Function to find tert-butyl hydrogens
def find_tert_butyl_hydrogens(atomic_info, bonds):
    """
    Finds all hydrogen atoms belonging to tert-butyl groups based on bonding information.
    
    :param atomic_info: Dictionary {atom_id: [symbol, coordinates]}.
    :param bonds: List of bonded atom pairs (0-based indices).
    :return: List of 0-based indices of tert-butyl hydrogen atoms.
    """
    tert_butyl_hydrogens = []
    bond_graph = nx.Graph(bonds)
    
    for atom_id, (symbol, _) in atomic_info.items():
        if symbol == 'C':  # Check if the atom is a carbon
            neighbors = list(bond_graph.neighbors(atom_id))
            
            methyl_count = 0
            temp_H = []
            for neighbor in neighbors:
                second_neighbors = list(bond_graph.neighbors(neighbor))
            
                # Separate neighbors into carbons and hydrogens
                carbons = [n for n in second_neighbors if atomic_info[n][0] == 'C']
                hydrogens = [n for n in second_neighbors if atomic_info[n][0] == 'H']
                temp_H.extend(hydrogens)
                
                if len(carbons) == 1 and len(hydrogens) == 3:
                    methyl_count += 1
            
            if methyl_count == 3:
                tert_butyl_hydrogens.extend(temp_H)

    return tert_butyl_hydrogens

# Example Usage:
pao_file = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\PAOr\pre\PAOr_pre.xyz"  # Provide the actual file path

ao_file = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\A0003\pre\A0003_pre.xyz"    # Provide the actual file path

pao_atom_idx = 51       # Change as per your system
ao_atom_idx = 9# 11        # Change as per your system
# pao_r_atom = 52  # User-defined reactive atom for PAO
# ao_r_atom = 40  # User-defined reactive atom for AO

steric_hindrance, pao_neighbors, ao_neighbors = calculate_gamma(
    pao_file, ao_file, pao_atom_idx, ao_atom_idx)
print(f"Steric Hindrance Factor (Gamma): {steric_hindrance}")
print(f"PAO Neighboring Atom Indices: {len(pao_neighbors)}--{pao_neighbors}")
print(f"AO Neighboring Atom Indices: {len(ao_neighbors)}--{ao_neighbors}")


#%%
atomic_info = read_xyz(ao_file)
bonds = calculate_bonds(atomic_info)
coordinates = np.array([info[1] for info in atomic_info.values()]) 
atomic_symbols = [info[0] for info in atomic_info.values()]
    
##
G = nx.Graph()
G.add_edges_from(bonds)
fig,ax=plt.subplots(dpi=1200)
types = [{'H':1, 'C':2, 'O':3}[x] for x in atomic_symbols]
atom_types = dict(zip(range(0,len(coordinates)), types))
man.draw_rdkit2D_from_graph(G, atom_types, "HCO", ax=ax,
                            atom_label=False, highlight_atoms=ao_neighbors)
ax.set_title('A0003')