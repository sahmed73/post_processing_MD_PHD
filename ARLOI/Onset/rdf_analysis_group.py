# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 28 13:08:35 2024
"""

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(dpi=350)

AOs = ['A0001','A0002']
for AO in AOs:
    # Load the trajectory and topology files
    topo_file=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\{}\Equilibration\Sim-1\H_removed_eq_1.data".format(AO)
    trej_file=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\{}\Production\Sim-1\prod.nvt.unwrapped.lammpstrj".format(AO)
    u = mda.Universe(topo_file, trej_file, format='LAMMPSDUMP')
    
    # Select donor hydrogens in antioxidants
    donor_hydrogens = u.select_atoms('index 12:1760')[::92]
    
    # Select reactive sites in radicals (e.g., oxygen atoms)
    reactive_sites = u.select_atoms('index 1873:2601')[::52]
    
    
    if AO == 'A0002':
        select = []
        for i in range(15):
            # select.append(1885 + i * 64)
            select.append(1887 + i * 64)
            # select.append(1893 + i * 64)
            select.append(1895 + i * 64)

        select_str = ' '.join(map(str, select))
        donor_hydrogens = u.select_atoms(f'index {select_str}')
        
    
    # Calculate RDF between donor hydrogens and reactive sites
    rdf = InterRDF(donor_hydrogens, reactive_sites, nbins=75, range=(0.0, 30.0))
    rdf.run()
    
    ax.plot(rdf.bins, rdf.rdf, label=AO)
#%%
   
ax.set_xlabel('Distance (Å)')
ax.set_ylabel('g(r)')
ax.set_title('RDF of O-H ')
ax.legend()
# plt.yscale('log')
