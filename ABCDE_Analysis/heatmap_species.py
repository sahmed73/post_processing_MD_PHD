# -*- coding: utf-8 -*-
"""
Created on Sat May 20 22:32:31 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import sys

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\Production\1350\Sim-1'

filename    = 'bonds.reaxc'
bondfile    = directory+'\\'+filename

index = directory.find('ABCDE')
molecule = directory[index+6:index+7]
print('Molecule:',molecule)

nat = bfp.get_neighbours(bondfile)
neighbours,atomtypes = nat
atomsymbols = ['H','C','O']
#%%-------------------
import re
lignins = {'A':'H12C10O3',
           'B':'H34C26O4',
           'C':'H44C29O2',
           'D':'H30C19O3',
           'E':'H14C11O4',
           'e':'H62C30'}
pattern = r'\b\d{3,4}\b'
matches = re.findall(pattern, directory)
print(matches)
final_temp = int(matches[0])
print('final temperature: ',final_temp)
skipts = (final_temp-300)/4
print('skip',skipts)
topspec = [lignins[molecule],'O2']
print('Molecue',lignins[molecule])
savedir = '..\\python_outputs\\heatmaps'
title   = directory[directory.find('LAMMPS'):]+'\n\nMolecule: {}'.format(molecule)
nspecies = 20
exclude = [lignins[molecule],'O2']

order = []

bfp.plot_species_heatmap(atomtypes,atomsymbols,neighbours,title=title,
                            savedir=savedir, nspecies=nspecies,
                            topspec=topspec, skipts=skipts,
                            exclude=exclude,log=True,figsize=(10,8)
                            )