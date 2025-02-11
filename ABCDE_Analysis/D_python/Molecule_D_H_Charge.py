# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 13:18:47 2023

@author: Shihab
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
atominfo = bfp.parsebondfile(bondfile,ALL=True)
#%%--
neighbours = atominfo['neighbours']
bondorders = atominfo['bondorders']
atomtypes  = atominfo['atypes']
charge     = atominfo['charge']
nlp        = atominfo['nlp']

steps = list(neighbours.keys())
it = 1200 # initial temperature
skipts = (it-300)/4
bondlist = []

firststepneigh = neighbours[steps[0]]
# 1 H
# 2 C
# 3 O
bonds = {}
flag = 0
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
                if flag < 1:
                    H = child
                    flag+=1
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
        
        # O-C Bonds
        if len(children)==2 and atomtypes[children[0]]==2 and atomtypes[children[1]]==2:
            for child in children:
                bond      = (parent,child)
                schildren = firststepneigh[child]
                sctypes   = [atomtypes[x] for x in schildren]
                
                # O-C-lower
                if sctypes.count(3) == 1:
                    key  = 'O-C-lower'
                    if key not in bonds.keys():
                        bonds[key] = [bond]
                    else:
                        bonds[key].append(bond)
                
                # O-C-upper
                elif sctypes.count(3) == 2:
                    key  = 'O-C-upper'
                    if key not in bonds.keys():
                        bonds[key] = [bond]
                    else:
                        bonds[key].append(bond)
                # 
                else:
                    print('Dumb')
bond_type = ['O-H','C-C','O-C-lower']
marker = ['s','^','x','o','>','<']
color  = 'bgrkc'

charge_x = []
for step,neigh in neighbours.items():
    for parent, children in neigh.items():
        if parent == H:
            charge_x.append(charge[step][parent])
ps = bfp.step2picosecond(steps, 0.25)

plt.plot(ps,charge_x)

###-plot-####
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule D'
plt.title(title)
#plt.title('Number of bond vs time {}'.format(bond_type))
plt.legend()
plt.xlabel('Time (ps)')
plt.ylabel('Charge')
plt.ylim(0,1)
plt.savefig('..\\python_outputs\\bondplot\\bondplot-'+ne.randstr(), dpi=300,bbox_inches='tight')
