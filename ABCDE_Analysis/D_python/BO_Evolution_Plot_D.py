# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 13:18:47 2023

@author: Shihab
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne
import pandas as pd
import numpy as np

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
atominfo = bfp.get_neighbours(bondfile,bo=True)
#%%--
neighbours, atomtypes, bondorders = atominfo
steps = list(neighbours.keys())
bondlist = []

firststepneigh = neighbours[steps[0]]
# 1 H
# 2 C
# 3 O
bonds = {}
ethyl_oxy=[]
for parent,children in firststepneigh.items():    
    if atomtypes[parent] == 3:
        # C=O bond
        if atomtypes[parent]==3 and len(children)==1 and atomtypes[children[0]]==2:
            child = children[0]
            key = 'C=O'
            if key not in bonds.keys():
                bonds[key] = [(parent,child)]
            else:
                bonds[key].append((parent,child))
                
        # O-H bonds
        for child in children:
            if atomtypes[child]==1:
                key = 'O-H'
                if key not in bonds.keys():
                    bonds[key] = [(parent,child)]
                else:
                    bonds[key].append((parent,child))
            
        # C-C Bonds    
        if len(children)==1:
            key = 'C-C'
            child1 = children[0]
            for ch in firststepneigh[child]:
                if atomtypes[ch]==2:
                    bond = (child1,ch)
                    if key not in bonds.keys():
                        bonds[key] = [bond]
                    else:
                        bonds[key].append(bond)
        
        # C-O Bonds
        if len(children)==2 and atomtypes[children[0]]==2 and atomtypes[children[1]]==2:
            for child in children:
                bond      = (parent,child)
                schildren = firststepneigh[child]
                sctypes   = [atomtypes[x] for x in schildren]
                
                # C-O-lower
                if sctypes.count(3) == 1:
                    ethyl_oxy.append(child)
                    key  = 'C-O-lower'
                    if key not in bonds.keys():
                        bonds[key] = [bond]
                    else:
                        bonds[key].append(bond)
                
                # C-O-upper
                elif sctypes.count(3) == 2:
                    key  = 'C-O-upper'
                    if key not in bonds.keys():
                        bonds[key] = [bond]
                    else:
                        bonds[key].append(bond)
                # 
                else:
                    print('Dumb')
        # O=O Bonds
        if len(children)==1 and atomtypes[children[0]]==3:
            key = 'O=O'
            bond = (parent,children[0])
            if key not in bonds.keys():
                bonds[key] = [bond]
            else:
                bonds[key].append(bond)
                

key  = 'Ethyl-Oxygen Bond'
keyy = 'H Absorption'
species = 'H5C2'
bonds[key]  = []
bonds[keyy] = []
oxygens = []

# Ethyl-Oxygen Bond
for step,neigh in neighbours.items():
    for parent, children in neigh.items():
        if parent in ethyl_oxy:
            for child in children:
                bond = (parent,child)
                if child>2600 and bond not in bonds[key]:
                    bonds[key].append(bond)
                
                ############################
                schildren = neigh[child]
                for schield in schildren:
                    if atomtypes[schield] == 3 and schield not in oxygens:
                        oxygens.append(schield)
                ############################
# H Absorptionby Oxygen
print(len(oxygens))
for step,neigh in neighbours.items():
    for parent, children in neigh.items():
        if parent in oxygens:
            for child in children:
                bond = (parent,child)
                if atomtypes[child]==1 and bond not in bonds[keyy]:
                    bonds[keyy].append(bond)
print(len(bonds[keyy]))
#%%-
bond_type = ['O-H','C-O-lower','C-O-upper','C-C','C=O','O=O','Ethyl-Oxygen Bond','H Absorption']
marker = ['s','^','x','o','>','*']    
cutoff = {'O-H':1, 'C=O':1.5, 'C-C':1.4, 'C-O-upper':2, 'C-O-lower':1,'Ethyl-Oxygen Bond':1,'H Absorption':1}
tol    = {'O-H':0.9, 'C=O':0.9, 'C-C':0.9, 'C-O-upper':0.35, 'C-O-lower':0.9,'Ethyl-Oxygen Bond':0.7,'H Absorption':0.7}
color    = {'O-H':'b', 'C=O':'k', 'C-C':'c', 'C-O-upper':'m', 'C-O-lower':'r','Ethyl-Oxygen Bond':'g','H Absorption':'k'}
marker   = {'O-H':'s', 'C=O':'^', 'C-C':'x', 'C-O-upper':'o', 'C-O-lower':'>','Ethyl-Oxygen Bond':'o','H Absorption':'x'}
label   = {'O-H':'Bond-1', 'C=O':'Bond-5', 'C-C':'Bond-2', 'C-O-upper':'Bond-3', 'C-O-lower':'Bond-4','Ethyl-Oxygen Bond':'Ethyl-Oxygen Bonds','H Absorption':'H Absorption'}

plot_keys = ['O-H','C-C','C-O-upper','C-O-lower']
final_temp = 1200
skipts = (final_temp-300)/4
print('skip',skipts)

for key in plot_keys:
    b = bonds[key]
    title = '{} {}\nMolecule D\n'.format(label[key],key)
    savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution'
    steps = np.array(list(bondorders.keys()))
    timestep = 0.25
    bfp.bondorder_evolution(bondorders, b, ts=timestep,savedir=savedir,title=title,skipts=skipts,fontsize=12)
