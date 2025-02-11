# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Nov 19 01:19:23 2024
"""
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt

# File paths
trejdirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0001\Production\300-500K_TRR=1Kpps\Sim-1'
trej = trejdirr + r'\prod.nvt.unwrapped.lammpstrj'
topodirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0001\Equilibration\Sim-1'
topo = topodirr + r'\H_removed_eq_1.data'

# Load the Universe
u = mda.Universe(topo, trej, format='LAMMPSDUMP')

# Select atom groups
PAO_oxygen = u.select_atoms("type 3").select_atoms('resid 1:20')  # For example: oxygen atoms in PAO radicals
ANY_hydrogen = u.select_atoms("type 2")  # For example: hydroxyl hydrogen atoms in AO

# Calculate RDF
rdf = InterRDF(PAO_oxygen, ANY_hydrogen, nbins=100, range=(0.0, 10.0))  # Adjust the range as needed
rdf.run()

# Plot RDF
fig, ax = plt.subplots(dpi=300)
ax.plot(rdf.bins, rdf.rdf, label='RDF: Radical Oxygen - Any hydrogen', color='tab:blue', lw=1.5)
ax.set_xlabel('Distance (Å)', fontsize=12)
ax.set_ylabel('g(r)', fontsize=12)
ax.legend(fontsize=10)
fig.tight_layout()
plt.show()


