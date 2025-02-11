# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Feb 12 10:31:18 2024
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import numpy as np
import seaborn as sns
import random
from collections import Counter

baseoil = "PAO4"
simdir = ['Sim-1']#,'Sim-2','Sim-3']
data = []
bonddata = {}
for cutoff in np.array(range(30,80,5))/100:
    for sim in simdir:
        directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(baseoil,baseoil,sim)
        bondfilepath = directory+'\\bonds.reaxc'
        bonddata[sim] = bfp.parsebondfile(bondfilepath,cutoff=cutoff,
                                          mtypes=True,bo=True)

    sim = 'Sim-1'
    neighbors = bonddata[sim]['neighbours']
    atypes    = bonddata[sim]['atypes']
    mtypes    = bonddata[sim]['mtypes']
    bondorders= bonddata[sim]['bondorders']
    atomsymbols = 'HCO'
    
    over_coordination  = [0, 0, 0]
    under_coordination = [0, 0, 0]
    
    checklist = []
    for frame, (step, neigh) in enumerate(neighbors.items()):
        for parent, children in neigh.items():
            nchild = len(children)
            check = sorted([parent]+children)
            if check in checklist:
                continue
            checklist.append(check)
            
            
            if atypes[parent]==1:
                if nchild>1: over_coordination[0]+=1
                elif nchild<1: under_coordination[0]+=1
                
            elif atypes[parent]==2:
                if nchild>4: over_coordination[1]+=1
                elif nchild<4: under_coordination[1]+=1
                
            elif atypes[parent]==3:
                if nchild>2: over_coordination[2]+=1
                elif nchild<2: under_coordination[2]+=1
        
        # counting fragments
        if frame==325:
            N = 0
            molecules = bfp.get_molecules(neigh)
            for molecule in molecules:
                species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
                if species=='H62C30':
                    N+=1
                elif species!='O2':
                    print(species)
            n_frag = 25-N
    
    data.append(np.array([*over_coordination, *under_coordination, n_frag]))
    
    print(f'cutoff: {cutoff}, over_coordination: {over_coordination}, under_coordination: {under_coordination}\nNumber of fragments: {n_frag}')
#%%
npdata  = np.array(data)
cutoffs = np.array(range(30,80,5))/100

plt.plot(cutoffs,npdata[:,0],label='Hydrogen Over Coordination')
plt.plot(cutoffs,npdata[:,1],label='Carbon Over Coordination')
plt.plot(cutoffs,npdata[:,2],label='Oxygen Over Coordination')
plt.plot(cutoffs,npdata[:,3:6].sum(axis=1),label='Total Under Coordination')
plt.xlabel('Bond Order Cutoff',fontsize=15)
plt.ylabel('Number of unphysical bonds',fontsize=15)
plt.legend()
# plt.minorticks_on()
plt.tick_params(axis='both',which='both',labelsize='large')
plt.savefig(r'..\..\..\python_outputs\figures\n_unphysical.png', dpi=300, bbox_inches='tight')