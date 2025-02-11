# -*- coding: utf-8 -*-
"""
Created on Wed May 10 00:57:08 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import magnolia.needless_essential as ne
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import pickle

timestep = 0.25
print('Assigned Timestep: ',timestep)

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_300_O2\Production\Sim-1'

filename    = 'bonds.reaxc'
bondfile    = directory+'\\'+filename

# check if the bond pickle file exists
picklefile     = directory+ '\\' + 'pickle_bond.pkl'
if os.path.exists(picklefile):
    print('Bondfile is loading from Pickle...')
    load_timestep, output_get_neighbours = pickle.load(open(picklefile, 'rb'))
    if load_timestep!=timestep:
        sys.exit('Error: Assigned timestep does not match with loaded the timestep\nAssigned Timestep: {} and Loaded Timestep: {}'.format(timestep,load_timestep))
else:
    print('Calling get_neighbours function....')
    output_get_neighbours = bfp.get_neighbours(bondfile,atypes=True,bo=True)
    with open(picklefile,'wb') as pf:
        pickle.dump((timestep,output_get_neighbours), pf)

neighbours,atomtypes,bondorders = output_get_neighbours
atom_symbols = ['H','C','O']

#%%###############Specific for AO1#############

C_O_bond_1 = []
C_O_bond_2 = []
C_O_bond_3 = []
C_O_bond_4 = []

first_index = 0
for atom,neighbour in neighbours[first_index].items():
    a_type = atomtypes[atom]
    
    #neigh_type = list(map(lambda x:atomtypes[x],neighbour))
    
    if a_type==3:
        for neigh in neighbour:
            neigh_type = atomtypes[neigh]
            if neigh_type==2:
                second_neighbour = neighbours[first_index][neigh]
                second_neighbour_type = list(map(lambda x:atomtypes[x],second_neighbour))
                if len(second_neighbour)==4:
                    if second_neighbour_type.count(1)==3:
                        C_O_bond_3.append((atom,neigh))
                    else:
                        C_O_bond_4.append((atom,neigh))
                

print(len(C_O_bond_1))
print(len(C_O_bond_2))
print(len(C_O_bond_3))
print(len(C_O_bond_4))

print(C_O_bond_4)
    
    
    
    
    
    
    