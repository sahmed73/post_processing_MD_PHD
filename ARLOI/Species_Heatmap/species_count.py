# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Oct  4 11:57:07 2024
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt

### Directory ###
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Equilibration\Sim-1"
filename  = "\\bonds.out"
bondfile = directory+filename
timestep  = 0.25
species= bfp.get_species_count(bondfile, 'HCO', timestep=0.25,
                               time_as_index=True)
#%%
temp=species.index*4+300
# species.index=temp
plt.plot(species,label=bfp.make_molecular_formula_latex(species.columns,sort=True))
plt.legend()
print(species)