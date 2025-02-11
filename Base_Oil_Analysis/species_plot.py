# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Oct 12 19:40:18 2023
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd


### Directory ###
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1100K"
filename  = "\\bonds.reaxc"
bondfilepath = directory+filename

### Parsing Bondfile ###
bonddata = bfp.parsebondfile(bondfilepath)

#%%## Geting neighbour lists
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
asyms      = ['H','C','O']
steps      = np.array(list(neighbours.keys()))


whole   = "H62C30"
seek = ["H2O","O2"]
df = {}

## Looping through steps and molecules
for step,neigh in neighbours.items():
    molecules = bfp.get_molecules(neigh)
    df[step]={}
    for x in seek: df[step][x] = 0 
    
    for molecule in molecules:
        species = bfp.get_molecular_formula(molecule, atypes, asyms)
        if species in seek:
            df[step][species]+=1
#%%
df0 = pd.DataFrame(df)
df1 = df0.T
df1.index = df1.index*0.25/1000
df1.index.name = 'time'
speciesCount = df1.reset_index()
print(speciesCount)



sns.set_style("whitegrid")
## ploting
latex = bfp.make_molecular_formula_latex(seek)
for i,species in enumerate(seek):
    sns.lineplot(x='time', y=species, data=speciesCount, label=latex[i])

