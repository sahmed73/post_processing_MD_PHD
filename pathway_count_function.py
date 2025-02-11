# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 01:40:48 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1000K'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
atominfo = bfp.get_neighbours(bondfile,bo=True)
atominfo += (atomsymbols,)
#%%----------------------------------------------
neighbours,atomtypes,bondorders,atomsymbols = atominfo

species = 'CO2'

seeklist,seekstep = bfp.species_to_molecule(atominfo, species) 
#%%
# import matplotlib.pyplot as plt
# seeklist,seekstep = result        
# for i,s in enumerate(seeklist):
#     f = bfp.get_molecular_formula(s, atomtypes, atomsymbols)
#     print(seekstep[i])
#     print(f,s)
#     print()
# print(len(seeklist))
# # ps = bfp.step2picosecond(seekstep, 0.25)
# # plt.scatter(ps,[0]*len(ps),s=0.07)
# # plt.yticks([])