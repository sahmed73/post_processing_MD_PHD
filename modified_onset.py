# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 13:32:04 2023

@author: Shihab
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\Onset\Sim-1'
filename  = '\\bonds.reaxc'
bondfilepath = directory+filename

bonddata = bfp.parsebondfile(bondfilepath, bo=True)
#%%
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
bondorders = bonddata['bondorders']
steps      = neighbours.keys()
atomsymbols= ['H','C','O']
timestep   = 0.25

NC = 30 # number of carbon
NH = 62 # number of hydrogen
NO = 0  # number of oxygen
seekMW = NC*12.0107 + NO*15.9994 + NH*1.00794# molecular weight w/o H: NH*1.00794
ramp = 4
it   = 300

print('seekMW=',seekMW)
x,y = [],[]
for step,neigh in neighbours.items():
    molecules = bfp.get_molecules(neigh)
    count = 0
    for molecule in molecules:
        species   = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
        MW        = bfp.compute_molecular_weight(species)
        if seekMW == MW:
            print(it+ramp*step*timestep/1000,species,len(molecule))
            count+=1
    x.append(it+ramp*step*timestep/1000)
    y.append(count)


#%%
plt.scatter(x, y,color='black')
plt.ylim(-1,26)