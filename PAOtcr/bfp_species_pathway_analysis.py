# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Apr 18 02:15:12 2024
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt
mplt.custom_plot_features()

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3A_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K'
filename = '\\bonds.reaxc'
bondfile = dirr+filename
bonddata = bfp.parsebondfile(bondfile,bo=True,mtypes=True)
#%% Unique Molecule
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']

timestep  = 0.25
steps     = np.array(list(neighbors.keys()))
ps        = steps*timestep/1000

radicals = {'PAOtcr':'H61C30O', 'A':'H30C19O3'}

unique_molecules = []
for step, neigh in neighbors.items():
    molecules = bfp.get_molecules(neigh)
    for molecule in molecules:
        if molecule not in unique_molecules:
            unique_molecules.append(molecule)
#%% Tracking each unique molecule

for i, unique_mol in enumerate(unique_molecules):
    reactants = []
    for step, neigh in neighbors.items():
        molecules = bfp.get_molecules(neigh)
        for molecule in molecules:
            common = unique_mol & molecule
            if common==unique_mol:
                break
            elif common:
                species = bfp.get_molecular_formula(molecule, atypes, 'HCO')
                if species not in reactants:
                    reactants.append(species)
    
    print(i+1,bfp.get_molecular_formula(unique_mol, atypes, 'HCO'))
    for reac in reactants:
        print(reac)
    print('='*20)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    