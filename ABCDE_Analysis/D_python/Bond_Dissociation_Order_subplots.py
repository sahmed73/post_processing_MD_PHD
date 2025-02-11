# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 23:21:35 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne
import pandas as pd
import numpy as np

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
atominfo = bfp.parsebondfile(bondfile,bo=True)
#%%--
neighbours = atominfo['neighbours']
atomtypes  = atominfo['atypes']
bondorders = atominfo['bondorders']

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
for step,neigh in neighbours.items():
    for parent, children in neigh.items():
        if parent in oxygens:
            for child in children:
                bond = (parent,child)
                if atomtypes[child]==1 and bond not in bonds[keyy]:
                    bonds[keyy].append(bond)
#%%-
import random
bond_type = ['O-H','C-O-lower','C-O-upper','C-C','C=O','O=O','Ethyl-Oxygen Bond','H Absorption']
marker = ['s','^','x','o','>','*']    
cutoff = {'O-H':1, 'C=O':1.5, 'C-C':1.4, 'C-O-upper':2, 'C-O-lower':1,'Ethyl-Oxygen Bond':1,'H Absorption':1}
tol    = {'O-H':0.9, 'C=O':0.9, 'C-C':0.9, 'C-O-upper':0.35, 'C-O-lower':0.9,'Ethyl-Oxygen Bond':0.7,'H Absorption':0.7}
color    = {'O-H':'b', 'C=O':'k', 'C-C':'c', 'C-O-upper':'m', 'C-O-lower':'r','Ethyl-Oxygen Bond':'g','H Absorption':'k'}
marker   = {'O-H':'s', 'C=O':'^', 'C-C':'x', 'C-O-upper':'o', 'C-O-lower':'>','Ethyl-Oxygen Bond':'o','H Absorption':'x'}
label   = {'O-H':'Bond-1', 'C=O':'Bond-5', 'C-C':'Bond-2', 'C-O-upper':'Bond-3', 'C-O-lower':'Bond-4','Ethyl-Oxygen Bond':'Ethyl-Oxygen Bonds','H Absorption':'H Absorption'}

plot_keys = ['C-C']
final_temp = 1200
skipts = (final_temp-300)/4
print('skip',skipts)
ps = bfp.step2picosecond(steps, 0.25)

for key in plot_keys:
    b = bonds[key]
    savedir = '1000K_'
    steps = np.array(list(bondorders.keys()))
    timestep = 0.25
    # bfp.bondorder_evolution(bondorders, b, ts=timestep,savedir=savedir,title=title,skipts=skipts,figsize=(8,10))
    
    result = {}
    beads  = dict(zip(range(1,51),[0]*50))
    cutoff = 0.7
    flag   = dict(zip(range(1,51),[0]*50))
    for step, bondorder in bondorders.items():
        pico = step*0.25/1000-skipts
        for i,atompair in enumerate(b):
            atom1,atom2 =atompair
            bo = bondorder[atom1].get(atom2,1)
            if bo<cutoff and flag[i+1]<3:
                beads[i+1]=pico
                flag[i+1] +=1
            if i+1 not in result.keys():
                result[i+1]=[bo]
            else:
                result[i+1].append(bo)

    print(len(beads),beads)
    ## shifting
    cut_index = ps.index(skipts)
    print(cut_index)
    dif = {}
    for k in result.keys():
        r = result[k]
        result[k] = np.array(r[cut_index:])
        dif[k]    = np.diff(r[cut_index-1:])
    npps = np.array(ps)
    ps = npps[npps>=skipts]-skipts
    
    #plotting
    plt.style.use('classic')
    fig, ax = plt.subplots(nrows=10, ncols=5)
    fig.set_size_inches(25.5, 28.5)
    j = 1
    for row in ax:
        for col in row:
            col.scatter(ps, result[j],s=0.6,label='ID: {}'.format(j),color='grey')
            col.plot([beads[j]]*10,np.linspace(cutoff,cutoff+0.5,10),color='red',linewidth=3)
            col.set_ylim(0.5,2)
            col.set_xlim(-5,1000)
            # col.set_xticks([0,250,500,750,1000])
            col.legend()
            j+=1
    plt.suptitle('{} {}\nMolecule D\n'.format(label[key],key),fontsize=35)
    fig.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\Bond Dissociation Order\D_Bond_Subplots\subplots-'+str(random.randint(0,10000000)), dpi=300)
    plt.show()

