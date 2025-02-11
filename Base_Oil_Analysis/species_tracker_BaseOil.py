# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Oct 28 23:03:17 2023
"""

################################
##
#   Species tracker
#   I want to make a tree of reaction
##
################################

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
bonddata = bfp.parsebondfile(bondfilepath,mtypes=True)
#%% Assigning Molecule Identification Number
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
atomsymbols= 'HCO'
mtypes     = bonddata['mtypes']

def get_label(MIN,search):
    label = next((k for k,v in MIN.items() if v==search), None)
    return label
    

molin_map      = {} # MOLecule Identification Number. key=molin,value=molecule
frame_map      = {}

molin         = 1

## Assigning Molecule Identification Number
for step, neigh in neighbours.items():
    molecules = list(bfp.get_molecules(neigh))
    for molecule in molecules:
        if molecule not in molin_map.values():
            molin_map[molin]=molecule.copy()
            molin+=1       
    

#%%# Tracking its parent
frame = 0
lastMolecules = []
parent_map     = {}
for step, neigh in neighbours.items():
    molecules = list(bfp.get_molecules(neigh))
    for molecule in molecules:
        parent = set()
        molin = get_label(molin_map,molecule)
        species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
        if molin not in parent_map.keys():
            for eachLastMol in lastMolecules:
                if eachLastMol==molecule:
                    continue
                intersection = eachLastMol & molecule
                if intersection:
                    last_molin = get_label(molin_map, eachLastMol)
                    parent.add(last_molin)
            
            parent_map[molin]=parent
    lastMolecules=molecules.copy()
    frame+=1
#%%
pathways = {}
seek = 'H2CO2'

def DFS(G, start_node):
    visited = set()  # Use a set for constant time lookups
    stack = [start_node]
    pathway = []

    while stack:
        v = stack.pop()
        if v not in visited:
            pathway.append(v)
            visited.add(v)  # Add to visited
            for w in G[v]:  # Iterate over neighbors
                if w not in visited:
                    stack.append(w)

    return pathway

f = open('result.txt','w')
for molin, molecule in molin_map.items():
    species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
    if species==seek:
        print(molecule)
        pathways[molin] = DFS(parent_map,molin)
        
for key, value in pathways.items():
    species = bfp.get_molecular_formula(molin_map[key], atypes, atomsymbols)
    print(species,file=f)
    for v in value:
        molecule = molin_map[v]
        u = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
        print(v,u,end='-->',file=f)
    print('\n----------',file=f)
    print(file=f)
f.close()