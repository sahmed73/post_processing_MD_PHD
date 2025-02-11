# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 07:33:39 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import magnolia.speciesfile_parser as sfp

directory    = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\200_Less_TempRamp_1150K'
filename     = '\\bonds.reaxc'
bondfilepath = directory+filename
atomsymbols  = ['H','C','O']
bondinfo = bfp.parsebondfile(bondfilepath, bo=True)
#%%
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

#%%
import magnolia.needless_essential as ne
import magnolia.speciesfile_parser as sfp
import sys
marker = ['s','^','x','o','>','*']
color    = ['b','c','purple','r','g','grey']
cutoff = 1
tol    = 0.4
skipts = (1150-300)/4
markevery = 150

species = sfp.get_species_count(directory+'\\species.out')
# species.index = species.index*0.25/1000
# species = species.loc[skipts:,:]
# species.index = species.index-skipts
print(species.loc['H2O'].values.size)
bonds['OH-H']=species.loc['H2O'].values[-4000:]
print(bonds['OH-H'])

subtitle = ["O–H bond","C–O bond","C–C bond","C=C bond"]
bondname = {'O-H':'1B: {} dissociation'.format(subtitle[0]),
            '-OH':'2B: {} dissociation'.format(subtitle[1]),
            'link':'3B: {} dissociation'.format(subtitle[2]),
            'C=C':'4B: {} dissociation'.format(subtitle[3]),
            'OH-H': 'HO-H bond formation'
            }

for i,key in enumerate(bondname.keys()):
    if key in ['OH-H']:
        plt.plot(ps,bonds[key],label='HO-H bond formation',marker='>',
                 markevery=markevery,markersize=4,linewidth=0.8,
                 color=color[i])
        continue
    print(key,len(bonds[key]))
    if key in ['C=C']:
        cutoff = 1.5
        tol    = 0.4
    ps, nbonds = bfp.get_nbondsVStime(bondorders, atypes, bonds[key],
                                      skipts=skipts,cutoff=cutoff, tol=tol)
    if key=='C=C': 
        nbonds = nbonds/2
    plt.plot(ps,nbonds,label=bondname[key],marker=marker[i],
             markevery=markevery,markersize=4,linewidth=0.8,
             color=color[i])

###-plot-####
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule B\n'
plt.title(title)
#plt.title('Number of bond vs time {}'.format(bond_type))
plt.legend(bbox_to_anchor=(1.58,1.02))
plt.xlabel('Time (ps)',fontsize=15)
plt.ylabel('Number of bonds',fontsize=15)
plt.xlim(-1,1001)
# plt.ylim(-5,205)
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\bondplot'
# plt.grid('on')
plt.savefig(savedir+'//_B-'+ne.randstr(), dpi=300,bbox_inches='tight')