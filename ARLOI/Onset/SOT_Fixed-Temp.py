# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  1 04:54:42 2024

SOT = Scavenging Onset Temperature
"""

import numpy as np
import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt

dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Production\fixed-300K"
filename=r"\bonds.out"
bondfile=dirr+filename
atominfo=bfp.parsebondfile(bondfile,mtypes=True)
#%%

neighbours=atominfo["neighbours"]
atypes=atominfo["atypes"]
mtypes=atominfo["mtypes"]

H, C, O = 1, 2, 3 # atom types
Lmin = 30 # minimum chain length
Lmax = 30 # maximum chain length
timestep = 0.25
N_PAOr = 20 # number of radicals
N_AO = 15 # number of antioxidants

SR = {} # count of scavenged radicals
for step, neigh in neighbours.items():
    molecules=bfp.get_molecules(neigh)
    G = nx.Graph(neigh)
    time = step*timestep/1000 # in ps
    SR[time] = 0
    for molecule in molecules:
        g = G.subgraph(molecule)
        L = [atypes[x] for x in molecule].count(C) # only counting carbons
        if L<Lmin or L>Lmax: continue
        for parent in g.nodes():
            c_types  = sorted([atypes[x] for x in g[parent]])
            if atypes[parent]==O and mtypes[parent]<=N_PAOr and c_types==[H,C]:
                SR[time]+=1
#%%
fig, ax =plt.subplots(dpi=350)    
ax.plot(SR.keys(), SR.values())
ax.set_xlabel("Time (ps)")
ax.set_ylabel("Number of SR")