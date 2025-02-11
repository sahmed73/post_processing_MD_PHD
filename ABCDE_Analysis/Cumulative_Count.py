# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 01:34:58 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import random 

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K'
filename    = 'bonds.reaxc'
bondfilepath= directory+'\\'+filename

bonddata = bfp.parsebondfile(bondfilepath,pkl='yes')
#%%------
neighbours  = bonddata['neighbours']
atypes      = bonddata['atypes']
atomsymbols = ['H','C','O']

species     = ['H26C16O3']
skipts      = (1200-300)/4

ps, ccount = bfp.cumulative_nspecies(neighbours, atypes, atomsymbols,
                                 species,skipts=skipts)
#%%
molecule = directory[directory.find('ABCDE')+8]
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\cumulative_species_count'
savedir = savedir + '\\cumulative_count_{}_{}'.format(species,random.randint(999,1000000))
title = directory[directory.find('ABCDE'):]+'\n\nMolecule {}'.format(molecule)  
plt.plot(ps,ccount[0],label= bfp.make_molecular_formula_latex(species,sort=True))
plt.title(title)
plt.xlabel('Time (ps)')
plt.ylabel('Cumulative number of molecules')
plt.legend()
plt.savefig(savedir, dpi=400,bbox_inches='tight')
print('--')