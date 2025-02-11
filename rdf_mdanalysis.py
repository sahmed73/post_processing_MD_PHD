# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Jan 22 15:14:30 2025
"""

# this is a test program to understand the rdf using mdanalysis
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt

# Load the universe
topofile = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0003\DataFile\Bulk\100_PAOr_50_A0003_density=0.2.data'
trejfile = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0003\Reaction\Sim-6_TSD=1.90\reaction.unwrapped.lammpstrj'
u = mda.Universe(
    topofile, trejfile,
    format="LAMMPSDUMP",
    custom_cols=["id", "resid", "type", "x", "y", "z"]  # Column names in the dump file
)
#%%
# Select reference oxygen atoms (in PAO radicals)
ref_oxygens = u.select_atoms("resid 1:100 and type 8")

# Select hydroxyl hydrogen atoms (in antioxidants)
hydroxyl_hydrogens = u.select_atoms("resid 101:150 and type 7")  # Adjust the selection as needed

#%%

# Compute the RDF
rdf = InterRDF(ref_oxygens, hydroxyl_hydrogens, range=(0.0, 10.0), nbins=200)
rdf.run()

# Plot the RDF
fig, ax = plt.subplots(dpi=350)
ax.plot(rdf.bins, rdf.rdf, label="RDF of H around O")
plt.xlabel("Distance (Å)")
plt.ylabel("g(r)")
plt.title("Radial Distribution Function")
plt.legend()
plt.grid()
plt.show()