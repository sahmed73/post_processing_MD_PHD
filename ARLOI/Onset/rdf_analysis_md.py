# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 28 12:38:36 2024
"""

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
fig, ax = plt.subplots(dpi=350)

# Load the trajectory and topology files
topo_file=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Equilibration\Sim-1\H_removed_eq_1.data"
trej_file=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Production\Sim-1\prod.nvt.unwrapped.lammpstrj"
u = mda.Universe(topo_file, trej_file, format='LAMMPSDUMP')
#%%

# Select donor hydrogens in antioxidants
donor_hydrogens = u.select_atoms('index 12:1760')[::92]

# Select reactive sites in radicals (e.g., oxygen atoms)
reactive_sites = u.select_atoms('index 1873:2601')[::52]

# Calculate RDF between donor hydrogens and reactive sites
rdf = InterRDF(donor_hydrogens, reactive_sites, nbins=75, range=(0.0, 20.0))
rdf.run()

# Plotting the RDF
plt.figure(dpi=350)
plt.plot(rdf.bins, rdf.rdf, label='RDF (H-O)')
plt.xlabel('Distance (Å)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function (RDF)')
plt.legend()
plt.show()
