# -*- coding: utf-8 -*-
"""
Created on Sat May 20 22:32:31 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import sys

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Squalane\25_Squalane_200_O2_Soria\Production\1600\Sim-2'
filename    = 'bonds.reaxc'
bondfile    = directory+'\\'+filename

bondinfo = bfp.parsebondfile(bondfile)
atomsymbols = ['H','C','O']
#%%-------------------
neighbours = bondinfo['neighbours']
atypes    = bondinfo['atypes']

savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\heatmaps'

title   = directory[directory.find('LAMMPS'):]+'\n'
nspecies = 15

order   = ['H2O', 'H2CO', 'H2O2', 'H61C30', 'H60C30O', 'H4C2O', 'CO2', 'HO2', 'H6C3', 'H6C3O', 'H34C18', 'H62C30O2', 'H10C5O', 'H44C22O2']
exclude = []#['O2','H62C30']
skipts = (1600-300)/4

bfp.plot_species_heatmap(atypes,atomsymbols,neighbours,
                         title=title,savedir=savedir,nspecies=nspecies,
                         order=order,sliced=[175,1000],
                         log=False,skipts=skipts,fontsize=14,exclude=exclude)
