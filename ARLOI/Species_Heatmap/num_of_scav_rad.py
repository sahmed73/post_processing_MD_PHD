# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Oct 13 03:12:28 2024
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import numpy as np
from pymatgen.core import Composition
import matplotlib.pyplot as plt

def count_scavenger_radicals(atominfo, chemical_formula, c_tol=0, h_tol=0):
    """
    Count the number of scavenger radicals over time based on atomic information.
    
    Parameters:
        atominfo (dict): Dictionary containing atomic information (neighbors, atypes, mtypes) for each time step.
        chemical_formula (str): The chemical formula of the molecule (e.g., "C6H12O6").
        c_tol (int, optional): Allowed difference between the initial and current carbon count. Default is 0.
        h_tol (int, optional): Allowed difference between the initial and current hydrogen count. Default is 0.
    
    Returns:
        dict: A time series of scavenger radical counts (keys=step, value=count).
    """
    
    neighbours = atominfo['neighbours']
    atypes  = atominfo['atypes']
    mtypes = atominfo['mtypes']
    
    composition = Composition(chemical_formula)
    n_carbon = composition['C']
    n_hydrogen = composition['H']
    
    n_scav_rad = {}
    for step, neigh in neighbours.items():
        graph = nx.Graph(neigh)
        molecules=bfp.get_molecules(neigh)
        n_scav_rad[step]=0
        for molecule in molecules:
            # current number of carbon and hydrogen in this molecule
            curr_n_carbon = [atypes[x] for x in molecule].count(2)
            curr_n_hydrogen = [atypes[x] for x in molecule].count(1)
            
            subgraph = graph.subgraph(molecule)
            for node in subgraph.nodes:
                if atypes[node]==3 and mtypes[node]<=20:
                    children=list(subgraph.neighbors(node))
                    children_type=sorted([atypes[x] for x in children])
                    if (children_type == [1, 2] 
                        and curr_n_carbon >= n_carbon - c_tol 
                        and curr_n_hydrogen >= n_hydrogen - h_tol):
                        
                        n_scav_rad[step]+=1
                    
    
    return n_scav_rad
                
#%%
AO = 'A0002'
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_{}\Production\Sim-1".format(AO)
filename  = "\\bonds.reaxc"
bondfile = directory+filename
timestep  = 0.25

atominfo = bfp.parsebondfile(bondfile, mtypes=True)
#%%

n_scav_rad=count_scavenger_radicals(atominfo, 'C30H62O', c_tol=0, h_tol=0)

##%% plotting

x = 300+4*np.array(list(n_scav_rad.keys()))/4000
y = np.array(list(n_scav_rad.values()))

plt.rcParams["font.size"]=15
fig, ax = plt.subplots(dpi=350)
ax.scatter(x,y)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Number of scavanged radical')

