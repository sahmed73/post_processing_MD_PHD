# -*- coding: utf-8 -*-
"""
Created on Thu May 25 01:51:34 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt

timestep = 0.25
print('Timestep: ',timestep)

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_400_O2\Production\TempRamp'
filename    = 'bonds.reaxc'
bondfile    = directory+'\\'+filename

nat = bfp.get_neighbours(bondfile)
neighbour,atomtypes = nat
atomsymbols = ['H','C','O']
#%%
species = 'H30C19O3'
ps,species_count=bfp.get_speciesVStemp(species,*nat,atomsymbols)
plt.plot(ps,species_count)
