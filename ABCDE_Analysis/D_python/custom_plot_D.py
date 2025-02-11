# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 02:01:01 2023

@author: arup2
"""

import matplotlib.pyplot as plt
import numpy as np
import magnolia.bondfile_parser as bfp
import sys
import random

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\Sim-2'
filename    = 'bonds.reaxc'
bondfile_path    = directory+'\\'+filename

bonddata = bfp.parsebondfile(bondfile_path, bo=True, mtypes=True)

#%%--
plt.style.use('fivethirtyeight')
plt.figure().set_figwidth(12)
def find_data_change_indices(data_list, change_value):
    change_indices = 0
    for i in range(len(data_list) - 1):
        if abs(data_list[i+1] - data_list[i]) >= change_value:
            # print(data_list[i+1] , data_list[i])
            # print('hi')
            change_indices = i+1
            break
    
    return change_indices

neighbours = bonddata['neighbours']
atomtypes  = bonddata['atypes']
mtypes     = bonddata['mtypes']
bondorders = bonddata['bondorders']

steps      = list(neighbours.keys())
picosecond = bfp.step2picosecond(steps, 0.25)

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

title = 'Bond Cleavage Order: Molecule D'


keys   = ['O-H','C-C','C-O-upper','C-O-lower']
labels= ['Bond-1','Bond-2','Bond-3','Bond-4']
colors = ['maroon','green','blue','black']
widths = [0.9,0.65,0.4,0.15]

first_skip = True

striks = []
for ii,key in enumerate(keys):
    bo = {}
    for i, bond in enumerate(bonds[key]):
        atom1, atom2 = bond
        bo[i+1] = []
        for step,ps in zip(steps,picosecond):
            bo[i+1].append(bondorders[step][atom1].get(atom2,0))
    # print(bo)       
    
    bits = []
    for bn, bo_value in bo.items():
        change_index = find_data_change_indices(bo_value,0.5)
        bits.append(picosecond[change_index])
    bits = [max(x-skipts,0) for x in bits]
    print(len(bits))
    
    #############-Sort-###############################################
    # if first_skip:
    #     firstbits = bits.copy()
    #     first_skip = False
    
    # if not first_skip:
    #     bits = sorted(bits,key = lambda x: firstbits[bits.index(x)])
    ##################################################################
           
    
    n = 100
    nmols = len(bits)
    tol   = 0.25
    y_top = 1000
    
    # blue line
    # for i in range(1,nmols+1):
    #     y = np.linspace(0,y_top,n)
    #     x = np.array([i]*n)
        
        # plt.plot(x,y,color='skyblue',alpha=0.9)
    
    # bits
    # for  i in range(1,nmols+1):
    #     x = np.linspace(i-tol, i+tol,n)
    #     y = np.array([bits[i-1]]*n)
    #     plt.plot(x,y,color=colors[ii])
    # # plt.bar(range(1,1+nmols),bits,color=colors[ii],width=widths[ii])
    # plt.bar(range(1,1+nmols),bits,facecolor='grey',width=0.6,alpha=0.1)
    # plt.bar(range(1,1+nmols),bits,facecolor='none',width=0.6,edgecolor='black',linewidth=0.6)
    # plt.bar(range(1,1+nmols),[1000]*len(bits),facecolor='grey',width=0.6,edgecolor='black',linewidth=0.6,alpha=0.1)
    
    striks.append(bits)
for i,color in enumerate(colors):
    plt.plot([],[],color=color,label=labels[i])

################################
ticksize = 7
################################

plt.title(title,fontsize=10)
plt.xticks(range(1,1+nmols),fontsize=ticksize)
plt.yticks(range(0,y_top+50,100),fontsize=ticksize)
plt.legend(fontsize=ticksize+2)

plt.text(7,940,'$3\\rightarrow 4 \\rightarrow 2 \\rightarrow 1: 22\%$          $3\\rightarrow 4 \\rightarrow 1 \\rightarrow 2: 14\%$          $3\\rightarrow 2 \\rightarrow 4 \\rightarrow 1: 10\%$          Others: $54\%$',color='k',fontsize=12)

savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\custom_BOE\CBOE-{}'.format(str(random.randint(999,1000000)))
plt.savefig(savedir,dpi=400, bbox_inches='tight')


#Analysis
for i in range(4):
    result = list(zip(striks[0],striks[1],striks[2],striks[3]))

cmpr = list(zip([1]*50,[2]*50,[3]*50,[4]*50))

sorted_result = []
single = []
for a,b in zip(result,cmpr):
    c = sorted(b,key=lambda x: a[b.index(x)])
    sorted_result.append(c)
    if c not in single:
        single.append(c)
    
sorted_ = sorted(sorted_result)
print(sorted_)
print(len(sorted_))
print(len(single))

for s in single:
    c = sorted_.count(s)*2
    print(s,'--',c,'%')
# set_    = set(sorted_)

for s_ in sorted_:
    s_.remove(2)
    s_.remove(3)
print(sorted_)
print(len(sorted_))
print('\n[1,4]--{}%'.format(sorted_.count([1,4])*2))
print('\n[4,1]--{}%'.format(sorted_.count([4,1])*2))

count = 0
for r in result:
    if r.count(0)==0:
        count+=1
print(count)












  