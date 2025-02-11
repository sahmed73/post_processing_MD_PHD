# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 03:45:33 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\200_Less_TempRamp_1150K'
filename  = 'bonds.reaxc'
bondfilepath = directory+'\\'+filename

bondinfo = bfp.parsebondfile(bondfilepath,bo=True)
#%%------
import sys
neighbours = bondinfo['neighbours']
bondorders = bondinfo['bondorders']
atypes     = bondinfo['atypes']
steps      = list(neighbours.keys())
firststepneigh = neighbours[steps[0]]

bonds = {}
for parent,children in firststepneigh.items():
    #O-H bond
    key = 'O-H'
    if atypes[parent]==3:
        for child in children:
            if atypes[child]==1:
                bond = (parent,child)
                if key in bonds:
                    bonds[key].append(bond)
                else:
                    bonds[key]=[]
    
    #tert-butyl
    key = 'tert-butyl'
    if atypes[parent]==2:
        ff = [atypes[x] for x in children].count(2)
        if ff==4:
            for child in children:
                fff = [atypes[x] for x in firststepneigh[child]].count(2)
                if fff==3:
                    bond = (parent,child)
                    if key in bonds:
                        bonds[key].append(bond)
                    else:
                        bonds[key]=[]
    
    #link
    key = 'link'
    if atypes[parent]==2:
        ff = [atypes[x] for x in children].count(2)
        if ff==3:
            for child in children:
                fff = [atypes[x] for x in firststepneigh[child]].count(2)
                if fff==3:
                    bond = (parent,child)
                    if key in bonds:
                        bonds[key].append(bond)
                    else:
                        bonds[key]=[]
    #C-H
    key = 'C-H'
    if atypes[parent]==1:
        for child in children:
            if atypes[child]==2:
                bond = (parent,child)
                if key in bonds:
                    bonds[key].append(bond)
                else:
                    bonds[key]=[]
                    
    #C=C
    key = 'C=C'
    if atypes[parent]==2:
        Noxy = [atypes[x] for x in children].count(3)
        Nhyd = [atypes[x] for x in children].count(1)
        if Noxy==1 and Nhyd==2:
            for child in children:
                if atypes[child]==2:
                    for schild in firststepneigh[child]:
                        if atypes[schild]==2 and schild!=parent:
                            bond = (schild,child)
                            if key in bonds:
                                bonds[key].append(bond)
                            else:
                                bonds[key]=[]
    
    # -OH
    key = '-OH'
    if atypes[parent]==3:
        for child in children:
            if atypes[child]==2:
                bond = (parent,child)
                if key in bonds:
                    bonds[key].append(bond)
                else:
                    bonds[key]=[]

for b in bonds:
    print(b,len(bonds[b]))
#%%    
bondname = {'O-H':'Bond-1 dissociation',
            '-OH':'Bond-2 dissociation',
            'link':'Bond-3 dissociation',
            'C=C':'Bond-4 dissociation',
            #'tert-butyl':'NoNeed'
            }

skipts = (1150-300)/4
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution'
molecule = directory[directory.find('ABCDE')+8]
print(molecule)

for key in bondname:
    title  = directory[directory.find('ABCDE'):]+'\nMolecule {}: {}--{}\n'.format(molecule,key,bondname[key])
    bfp.bondorder_evolution(bondorders, bonds[key],skipts=skipts,title=title,
                        savedir = savedir,fontsize=12)