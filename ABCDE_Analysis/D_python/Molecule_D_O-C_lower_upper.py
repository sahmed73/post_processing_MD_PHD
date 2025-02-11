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

marker = ['s','^','x','o','>','<']
color  = 'b g r k c'.split()
index = 0       
cutoff = [1,1.5,1.4,1,1]
tol    = [0.9,0.9,0.9,1,0.9]
print(bonds.keys())
for key, bondlist in bonds.items():
    markevery = 300
    if key == 'O-C-upper': 
        cutoff[index] = 0.2
        tol[index] = 1.4
    
        
        
        ps, nbonds = bfp.get_nbondsVStime(bondorders,atomtypes,
                                          bondlist,skipts=skipts,
                                          cutoff=cutoff[index], tol=tol[index])
        
        plt.plot(ps,nbonds,label='Bond-3 Single Bond',marker=marker[index],
                 markevery=markevery,markersize=4,linewidth=0.8,
                 color=color[index])
        
        cutoff[index] = 2 #1.8
        tol[index] = .35 #0.5
        color[index] = 'm'
        marker[index] = 's'
        
        ps, nbonds = bfp.get_nbondsVStime(bondorders,atomtypes,
                                          bondlist,skipts=skipts,
                                          cutoff=cutoff[index], tol=tol[index])
        
        plt.plot(ps,nbonds,label='Bond-3 Double Bond',marker=marker[index],
                 markevery=markevery,markersize=4,linewidth=0.5,
                 color=color[index])
    
    index+=1

###-plot-####
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule D'
plt.title(title)
#plt.title('Number of bond vs time {}'.format(bond_type))
plt.legend()
plt.xlabel('Time (ps)')
plt.ylabel('Number of bonds')
plt.savefig('..\\python_outputs\\bondplot\\bondplot-'+ne.randstr(), dpi=300,bbox_inches='tight')
