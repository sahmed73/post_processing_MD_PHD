# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 04:16:32 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import magnolia.needless_essential as ne


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\C\C_300_O2\Production\LessThan_200K_1200'
filename  = 'bonds.reaxc'
bondfilepath = directory+'\\'+filename

bondinfo = bfp.parsebondfile(bondfilepath,mtypes=True,bo=True)
#%%------
neighbours = bondinfo['neighbours']
atypes     = bondinfo['atypes']
mtypes     = bondinfo['mtypes']
bondorders = bondinfo['bondorders']
steps      = list(neighbours.keys())


bondname = {'OH':'Bond-1','C-C_link':'Bond-2'}

####-Bond Extracting-##############
firststepneigh = neighbours[steps[0]]
# 1 H
# 2 C
# 3 O
bonds = {}
for parent,children in firststepneigh.items():
    # O-H bonds    
    if atypes[parent] == 3:                
        for child in children:
            if atypes[child]==1:
                key = 'O-H'
                if key not in bonds.keys():
                    bonds[key] = [(parent,child)]
                else:
                    bonds[key].append((parent,child))
    
    # C-C_link bonds
    if atypes[parent] == 2: # Carbon
        H_count = [atypes[x] for x in children].count(1) # Hydrogen
        C_count = [atypes[x] for x in children].count(2) # Carbon
        if H_count == 2 and C_count == 2:
            for child in children:
                if atypes[child]==2:
                    key = 'C-C_link'
                    if key not in bonds.keys():
                        bonds[key] = [(parent,child)]
                    else:
                        bonds[key].append((parent,child))
    
    # Tertiary Butyl link bonds
    if atypes[parent] == 2: #Carbon
        key = 'tertiary'
        if [atypes[x] for x in children].count(2)==4:
            for child in children:
                if [atypes[x] for x in firststepneigh[child]].count(2)==3:
                    
                    if key not in bonds.keys():
                        bonds[key] = [(parent,child)]
                    else:
                        bonds[key].append((parent,child))


# adding missing O-H bond
missing = (256,300)
bonds['O-H'].append(missing)
print('added missing O-H bond',missing)

# delete extra C-C_link bonds
extra = [(3301,3302),(3301,3306),(3302,3301),(3306,3301)]
bond = bonds['C-C_link']
bond = [x for x in bond if x not in extra]
bonds['C-C_link'] = bond

key = 'O-H'

for b in bonds:
    print(b,len(bonds[b]))


skipts = (1200-300)/4
molecule = directory[directory.find('ABCDE')+8]
print(molecule)
title  = directory[directory.find('ABCDE'):]+'\nMolecule {}: {} bonds'.format(molecule,key)
bfp.bondorder_evolution(bondorders, bonds[key],skipts=skipts,title=title,sort=700)

for atom1,atom2 in bonds['O-H']:
    hydrogen = atom1 if atypes[atom1]==1 else atom2
    # print(atypes[hydrogen])