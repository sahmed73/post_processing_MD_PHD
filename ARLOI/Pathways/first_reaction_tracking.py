# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Sep 13 11:18:02 2024
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt
import time
import magnolia.needless_essential as ne
import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import magnolia.access_ucm_cluster as ucm


directory = "/mnt/borgstore/amartini/sahmed73/ARLOI-V1/4_ARLOI_24-09-07--22-34-03/20_PAO-OH_15_A0001/Production/Sim-1"
filename  = "/bonds.reaxc"
bondfile = ucm.local_copy(directory+filename)
timestep  = 0.25

atominfo = bfp.parsebondfile(bondfile, bo=True)

neighbours = atominfo['neighbours']
atomtypes  = atominfo['atypes']
bondorders = atominfo['bondorders']
#%%

species = set()
for step, neigh in neighbours.items():
    molecules = bfp.get_molecules(neigh)
    for molecule in molecules:
        spec=bfp.get_molecular_formula(molecule, atomtypes, 'HCO')
        species.add(spec)
    if(len(species)>2):
        break
print(len(species),species, step)