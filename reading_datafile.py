# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 25 03:05:55 2024
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.rms import rmsd

dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\TEST\T1\Geometry\bulk-05\DataFile'
topofile=dirr+'\\50_BHT_100_PAOr_10_Rep.data'

# Load the LAMMPS data file
u = mda.Universe(topofile)

# Group molecules by residue ID
molecules = [res for res in u.residues]

def get_bond_connectivity(molecule):
    """
    Extracts bond connectivity as a set of tuples (atom1_id, atom2_id).
    """
    return {(bond.atoms[0].id, bond.atoms[1].id) for bond in molecule.bonds}

distinct_structures = []
for mol1 in molecules:
    is_distinct = True
    bonds1 = get_bond_connectivity(mol1)
    for mol2 in distinct_structures:
        bonds2 = get_bond_connectivity(mol2)
        if bonds1 == bonds2:  # Compare bond connectivity
            is_distinct = False
            break
    if is_distinct:
        distinct_structures.append(mol1)

print(f"Number of distinct structures: {len(distinct_structures)}")

