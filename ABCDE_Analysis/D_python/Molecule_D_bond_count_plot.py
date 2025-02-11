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

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
atomsymbols = ['H','C','O']
atominfo = bfp.get_neighbours(bondfile,bo=True)
#%%--
neighbours, atomtypes, bondorders = atominfo
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
                

key  = 'Ethyl-Oxygen Bond'
keyy = 'H Absorption'
species = 'H5C2'
bonds[key]  = []
bonds[keyy] = []
oxygens = []

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
color    = {'O-H':'b', 'C=O':'grey', 'C-C':'c', 'C-O-upper':'m', 'C-O-lower':'r','Ethyl-Oxygen Bond':'g','H Absorption':'grey'}
marker   = {'O-H':'s', 'C=O':'^', 'C-C':'x', 'C-O-upper':'o', 'C-O-lower':'>','Ethyl-Oxygen Bond':'^','H Absorption':'v'}
label   = {'O-H':'1A: O–H bond dissociation', 'C=O':'Bond-5', 'C-C':'2A: C–C bond dissociation', 'C-O-upper':'3A: C–O bond single to double', 'C-O-lower':'4A: C–O bond dissociation','Ethyl-Oxygen Bond':'Ethyl-oxygen bond formation','H Absorption':'H absorption bond formation'}
markevery = 150

plot_keys = ['O-H','C-C','C-O-upper','C-O-lower',
             'Ethyl-Oxygen Bond','H Absorption']
bond_formation = ['Ethyl-Oxygen Bond','H Absorption']
subtitle = ["O–H bond","C–C bond","C–O bond","C–O bond"]
print(bonds.keys())
fig, ax1 = plt.subplots()
for key in plot_keys:
    bondlist = bonds[key]
    if key== bond_formation:
        ps, nbonds = bfp.get_nbondsVStime(bondorders, atomtypes,
                                          bondlist,skipts=skipts,
                                          cutoff=cutoff[key], tol=tol[key])
        ax2 = ax1.twinx()
        ax2.plot(ps,nbonds,label=label[key],marker=marker[key],
                 markevery=markevery,markersize=4,linewidth=0.8,
                 color=color[key])
        continue

    if key == 'C-O-upper':
        markevery = 247
    ps, nbonds = bfp.get_nbondsVStime(bondorders, atomtypes,
                                      bondlist,skipts=skipts,
                                      cutoff=cutoff[key], tol=tol[key])
    
    ax1.plot(ps,nbonds,label=label[key],marker=marker[key],
             markevery=markevery,markersize=4,linewidth=0.8,
             color=color[key])

###-plot-####
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule D'+"\n"
ax1.set_title(title,fontsize=10)
#plt.title('Number of bond vs time {}'.format(bond_type))
ax1.legend(bbox_to_anchor=(1.63,1.02))
ax1.set_xlabel('Time (ps)',fontsize=15)
ax1.set_ylabel('Number of bonds',fontsize=15)
ax1.set_ylim(0,55)
ax1.set_xlim(right=1001)
fig.savefig('..\\..\\python_outputs\\bondplot\\bondplot-'+ne.randstr(), dpi=300,bbox_inches='tight')
