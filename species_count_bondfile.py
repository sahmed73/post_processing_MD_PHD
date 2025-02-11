# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Apr 14 20:27:04 2024
"""
import magnolia.bondfile_parser as bfp

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\Pure_50B\Production\1050\Sim-1'
filename = '\\bonds.reaxc'
bondfile = dirr+filename

species = bfp.get_species_count(bondfile, atomsymbols=['H', 'C', 'O', 'Si'],
                                timestep=0.25, restart_time=True)
print(species)
#%%