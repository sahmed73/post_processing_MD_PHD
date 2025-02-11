# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  3 02:12:49 2023
"""

import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
import random
import numpy as np

### Directory ###
base_oil = 'PAO4'
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(base_oil,base_oil)
filename  = "\\bonds.reaxc"
bondfilepath = directory+filename

### Parsing Bondfile ###
bonddata = bfp.parsebondfile(bondfilepath)
#%% Assigning Molecule Identification Number
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
atomsymbols= 'HCO'
timestep = 0.25
seek = {}

f = open('result.txt','w')
frame = 0 
for step, neigh in neighbours.items():
    ps = step*timestep/1000
    molecules = bfp.get_molecules(neigh)
    for molecule in molecules:
        if seek==molecule:
            species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
            print("Frame {}, Time {} ps\n{}".format(frame,ps,species),file=f)
    frame+=1
f.close()