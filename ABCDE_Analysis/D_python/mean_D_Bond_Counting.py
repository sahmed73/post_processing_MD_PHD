# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 06:26:19 2023

@author: arup2
"""
import matplotlib.pyplot as plt
import numpy as np
import magnolia.bondfile_parser as bfp
import sys
import random

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\Sim-2'
filename    = 'bonds.reaxc'
bondfile_path    = directory+'\\'+filename

bonddata = bfp.parsebondfile(bondfile_path, bo=True)
#%%--
def find_data_change_indices(data_list, change_value):
    change_indices = 0
    for i in range(len(data_list) - 1):
        if abs(data_list[i+1] - data_list[i]) >= change_value:
            # print(data_list[i+1] , data_list[i])
            # print('hi')
            change_indices = i+1
            break
    
    return change_indices

neighbours = bonddata['neighbours']
atomtypes  = bonddata['atypes']
bondorders = bonddata['bondorders']

steps      = list(neighbours.keys())
picosecond = bfp.step2picosecond(steps, 0.25)

steps = list(neighbours.keys())
it = 1200 # initial temperature
skipts = (it-300)/4
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

keys   = ['O-H','C-C','C-O-upper','C-O-lower']
labels= ['Bond-1','Bond-2','Bond-3','Bond-4']

striks = []
for ii,key in enumerate(keys):
    bo = {}
    for i, bond in enumerate(bonds[key]):
        atom1, atom2 = bond
        bo[i+1] = []
        for step,ps in zip(steps,picosecond):
            bo[i+1].append(bondorders[step][atom1].get(atom2,0))
    # print(bo)       
    
    bits = []
    for bn, bo_value in bo.items():
        change_index = find_data_change_indices(bo_value,0.5) ####
        bits.append(picosecond[change_index])
    bits = [max(x-skipts,0) for x in bits]
    print(len(bits),end='  ')
    striks.append(bits)
print()

#Analysis
result = list(zip(striks[0],striks[1],striks[2],striks[3]))
cmpr = list(zip([1]*50,[2]*50,[3]*50,[4]*50))

sorted_result = []
single = []
for a,b in zip(result,cmpr):
    c = sorted(b,key=lambda x: a[b.index(x)])
    sorted_result.append(c)
    if c not in single:
        single.append(c)
    
sorted_ = sorted(sorted_result)

for s in single:
    c = sorted_.count(s)*2
    print(s,'--',c,'%')

for s_ in sorted_:
    s_.remove(2)
    s_.remove(3)
print(sorted_)
print(len(sorted_))
print('\n[1,4]--{}%'.format(sorted_.count([1,4])*2))
print('\n[4,1]--{}%'.format(sorted_.count([4,1])*2))

count = 0
for r in result:
    if r.count(0)==0:
        count+=1
print(count)












  