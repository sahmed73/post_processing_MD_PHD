# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Sep 15 21:05:03 2024
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import magnolia.access_ucm_cluster as ucm
import magnolia.molecular_analysis as man

dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0004\Production\Sim-1"
filename  = r"\bonds.reaxc"
# bondfile = ucm.local_copy(dirr+filename)
bondfile = dirr+filename
timestep  = 0.25

atominfo = bfp.parsebondfile(bondfile, bo=True)

neighbours = atominfo['neighbours']
atomtypes  = atominfo['atypes']
bondorders = atominfo['bondorders']
#%%
atomConnectivity = bfp.parsebondfile_asGraph(bondfile)
#%%
atomsymbols='HCO'
initial_number=35 # user_input
AOspecies = 'H29C19O3' # user_input

last_count=None
for step, graph in atomConnectivity.items():    
    dynamic_number=nx.number_connected_components(graph)
    if dynamic_number!=initial_number:
        if dynamic_number!=last_count:
            last_count=dynamic_number
        else:
            print(step, dynamic_number)
            
            molecules = bfp.get_molecules(graph)
            for molecule in molecules:
                # subgraph = graph.subgraph(molecule)
                species = bfp.get_molecular_formula(molecule, atomtypes, atomsymbols)
                print(species)
            break