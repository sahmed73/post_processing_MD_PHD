# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Oct 11 23:51:48 2023
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd


### Directory ###
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_with_Odd_Density\Production\1350\Sim-1"
filename  = "\\bonds.reaxc"
bondfilepath = directory+filename

### Parsing Bondfile ###
bonddata = bfp.parsebondfile(bondfilepath,mtypes=True)

#%%## Geting neighbour lists
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
mtypes     = bonddata['mtypes']
asyms      = ['H','C','O']
steps      = np.array(list(neighbours.keys()))

whole   = "H62C30"

## Looping through steps and molecules
for step,neigh in neighbours.items():
    molecules = bfp.get_molecules(neigh)
    for molecule in molecules:
        species = bfp.get_molecular_formula(molecule, atypes, asyms)
        if species == whole:
            whole_molecule = molecule
            break
    break

print(whole_molecule)
print(len(whole_molecule))
tracker = open('tracker.txt','w')
frame = -1-1050

commonAtomIDs = whole_molecule

for step,neigh in neighbours.items():
    molecules = bfp.get_molecules(neigh)
    print_flag = True
    frame+=1
    if frame==0:
        print(step*0.25/1000)
    for molecule in molecules:
        intersection = molecule & whole_molecule
        if intersection == whole_molecule:
            continue
        elif intersection:
            commonAtomIDs.update(molecule)
            if print_flag:
                print("Frame: {} ,  Steps: {}".format(frame,step),file=tracker)
                print_flag = False
            print(bfp.get_molecular_formula(molecule, atypes, asyms),file=tracker)
tracker.close()
print(commonAtomIDs)
print(len(commonAtomIDs))

mtyp = set()
for atom in commonAtomIDs:
    mtyp.add(mtypes[atom])
print(mtyp)
print(len(mtyp))