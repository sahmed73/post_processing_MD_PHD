# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 06:57:36 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import magnolia.needless_essential as ne
import sys
import random
import numpy as np

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K'
filename    = 'bonds.reaxc'
bondfile_path    = directory+'\\'+filename

neighbours,atomtypes,bondorders,mtypes = bfp.get_neighbours(bondfile_path, bo=True, mtypes=True)
steps = list(neighbours.keys())
#%%--
import matplotlib.pyplot as plt
bonds = []
atoms = []

firststepneigh = neighbours[steps[0]]

carbons = []

for parent,children in firststepneigh.items():
    if atomtypes[parent]==3 and len(children)==2 and atomtypes[children[0]]==2 and atomtypes[children[1]]==2:
        for child in children:
            schildren = firststepneigh[child]
            sctypes   = [atomtypes[x] for x in schildren]
            # O-C-lower
            if sctypes.count(3) == 1:
                carbons.append(child)
                
print(len(carbons))
for step,neigh in neighbours.items():
    for parent, children in neigh.items():
        if parent in carbons:
            #print(atomtypes[parent],mtypes[parent],end=' | ',sep=',')
            for child in children:
                bond = (parent,child)
                if child>2600 and bond not in bonds:
                    bonds.append(bond)

#%%        
print('number of bond: ',len(bonds))
figsize = (5,10)
final_temp = 1200
skipts = (final_temp-300)/4
#print('skip',skipts)

cutoff = 1
tol    = 0.7
label  = 'Ethyl-Oxygen Bond'
marker = 's'
markevery = 100
color = 'black'

ps, nbonds = bfp.get_nbondsVStime(bondorders, atomtypes, bonds,skipts=skipts,cutoff=cutoff, tol=tol)
        
plt.plot(ps,nbonds,label=label,marker=marker,
         markevery=markevery,markersize=4,linewidth=0.8,
         color=color)

###-plot-####
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule D'
plt.title(title)
#plt.title('Number of bond vs time {}'.format(bond_type))
plt.legend()
plt.xlabel('Time (ps)')
plt.ylabel('Number of bonds')
plt.savefig('..\\..\\python_outputs\\bondplot\\bondplot-'+ne.randstr(), dpi=300,bbox_inches='tight')