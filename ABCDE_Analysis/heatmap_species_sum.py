# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 02:33:45 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import sys
import os

d = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\C\C_300_O2\Production\1277K\Sim-'

neighbours_list = []
for sim in '123':
    print('Simulation-{}'.format(sim))
    directory = d
    if os.path.exists(d+sim):
        directory = d+sim
    filename    = 'bonds.reaxc'
    bondfile    = directory+'\\'+filename
    
    index = directory.find('ABCDE')
    molecule = directory[index+6:index+7]
    print('Molecule:',molecule)
    
    nat = bfp.get_neighbours(bondfile)
    neighbours,atomtypes = nat
    neighbours_list.append(neighbours)
    if not os.path.exists(d+sim):
        break
#%%-------------------
atomsymbols = ['H','C','O']
lignins = {'A':'H12C10O3',
           'B':'H34C26O4',
           'C':'H44C29O2',
           'D':'H30C19O3',
           'E':'H14C11O4'}
 
final_temp = 1000
skipts = (final_temp-300)/4
print('skip',skipts)
topspec = [lignins[molecule],'O2']
print('Molecue',lignins[molecule])
savedir = '..\\python_outputs\\heatmaps'
title   = directory[directory.find('LAMMPS'):]+'\n\nMolecule: {}'.format(molecule)
nspecies = 10

order = None
exclude = [lignins[molecule],'O2']

bfp.plot_species_heatmap(atomtypes,atomsymbols,*neighbours_list,title=title,
                         savedir=savedir,
                         nspecies=nspecies,order=order,
                         topspec=topspec, skipts=skipts, kind='sum',
                         exclude=exclude,log=False,
                         fontsize=15)