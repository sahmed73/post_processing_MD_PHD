# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 03:20:08 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\C\C_300_O2\Production\LessThan_000K_1400'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
bondinfo = bfp.parsebondfile(bondfile,bo=True)
#%%--
neighbours = bondinfo['neighbours']
atypes     = bondinfo['atypes']
bondorders = bondinfo['bondorders']
steps      = list(neighbours.keys())

####-Bond Extracting-##############
firststepneigh = neighbours[steps[0]]
# 1 H
# 2 C
# 3 O
bonds = {}

CH3_O = [] # to find bond-5
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
                        
                    ## -CH3 Bonds
                    for carbon in [parent,child]:
                        ff = firststepneigh[carbon]
                        if [atypes[x] for x in ff].count(2)==4:
                            ttchildren = firststepneigh[carbon]
                            for ttchild in ttchildren:
                                keyy = '-CH3'
                                fff = firststepneigh[ttchild]
                                if [atypes[x] for x in fff].count(1)==3:
                                    CH3_O.append(ttchild)
                                    if keyy not in bonds.keys():
                                        bonds[keyy] = [(carbon,ttchild)]
                                    else:
                                        bonds[keyy].append((carbon,ttchild))
### finding extra bonds #######
key = 'H2C=O'
bonds[key]=[]
for step, neigh in neighbours.items():
    for parent, children in neigh.items():
        if parent in CH3_O:
            for child in children:
                bond = (parent,child)
                if child>3750 and bond not in bonds[key]:
                    bonds[key].append(bond)
               
bondname = {'O-H':'Bond-1 dissociation','C-C_link':'Bond-2 dissociation','tertiary':'Bond-3 dissociation', 
            '-CH3':'Bond-4 dissociation','H2C=O':'C-O bond formation'}
    

# adding missing O-H bond
missing = (256,300)
bonds['O-H'].append(missing)
print('added missing O-H bond',missing)

# delete extra C-C_link bonds
extra = [(3301,3302),(3301,3306),(3302,3301),(3306,3301)]
bond = bonds['C-C_link']
bond = [x for x in bond if x not in extra]
bonds['C-C_link'] = bond

skipts = 70
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution'
molecule = directory[directory.find('ABCDE')+8]
print(molecule)
title  = directory[directory.find('ABCDE'):]+'\nMolecule {}: {} bonds'.format(molecule,key)
#%%-
bond_type = list(bondname.keys())
marker = ['s','^','x','o','>','*']    
cutoff = {'O-H':1, 'C=O':1.5, 'C-C':1.4, 'C-O-upper':2, 'C-O-lower':1,'Ethyl-Oxygen Bond':1,'H Absorption':1}
tol    = {'O-H':0.9, 'C=O':0.9, 'C-C':0.9, 'C-O-upper':0.35, 'C-O-lower':0.9,'Ethyl-Oxygen Bond':0.7,'H Absorption':0.7}

cutoff = 1
tol    = 0.3

plot_keys = bond_type
print(plot_keys)
print(bonds.keys())
for i,key in enumerate(plot_keys):
    if key in plot_keys:
        bondlist = bonds[key]
        markevery = 150
        if key == 'C-O-upper':
            markevery = 247
        ps, nbonds = bfp.get_nbondsVStime(bondorders, atypes, bondlist,skipts=skipts,cutoff=cutoff, tol=tol)
        if key=='H2C=O':
            savebond = nbonds
        plt.plot(ps,nbonds,label=bondname[key],marker=marker[i],
                 markevery=markevery,markersize=4,linewidth=0.8)

###-plot-####
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule C\n'
plt.title(title)
#plt.title('Number of bond vs time {}'.format(bond_type))
plt.legend(bbox_to_anchor=(1,1))
plt.xlabel('Time (ps)')
plt.ylabel('Number of bonds')
# plt.xlim(right=1150)
plt.savefig('..\\..\\python_outputs\\bondplot\\bondplot-'+ne.randstr(), dpi=300,bbox_inches='tight')
