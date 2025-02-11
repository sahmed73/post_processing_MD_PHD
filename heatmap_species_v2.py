# -*- coding: utf-8 -*-
"""
Created on Thu May 25 02:13:12 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import sys

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO2\50_AO2_300_O2\Production\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
bondpickle_path = directory+'\\bonds.pickle'

# first bondfile is directory to find pickle
result = bfp.loadpickle_or_execute(bondpickle_path,bfp.get_neighbours,bondfile,bo=True,mtypes=True)

#dropping bo and mtypes and adding atomsymbols
atominfo = tuple(result[:2])+(['H','C','O'],) 

neighbours,atomtypes,atomsymbols = atominfo

#%%
savedir = 'python_outputs\\heatmaps'
title   = directory[directory.find('LAMMPS'):]
nspecies = 20

species_pickle_path = directory+'\\species.pickle'
bfp.plot_species_heatmap(*atominfo,title=title,savedir=savedir,nspecies=nspecies,pickle=species_pickle_path)