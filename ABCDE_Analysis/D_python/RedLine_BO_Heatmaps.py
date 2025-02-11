# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 01:44:25 2023

@author: arup2
"""
one={1: 260.75, 2: 39.0, 3: 184.75, 4: 325.75, 5: 329.5, 6: 317.25, 7: 218.25, 8: 136.0, 9: 87.5, 10: 355.5, 11: 431.25, 12: 628.25, 13: 298.75, 14: 377.0, 15: 305.75, 16: 587.5, 17: 702.0, 18: 688.0, 19: 421.0, 20: 576.25, 21: 150.75, 22: 1000, 23: 315.5, 24: 112.5, 25: 190.75, 26: 477.75, 27: 386.25, 28: 803.25, 29: 119.0, 30: 163.25, 31: 608.25, 32: 282.5, 33: 408.0, 34: 338.5, 35: 568.5, 36: 443.5, 37: 355.0, 38: 185.0, 39: 323.75, 40: 112.5, 41: 39.5, 42: 571.0, 43: 227.75, 44: 426.75, 45: 33.25, 46: 334.0, 47: 24.75, 48: 242.75, 49: 376.75, 50: 624.75}

two={1: 616.0, 2: 1000, 3: 713.75, 4: 65.25, 5: 1000, 6: 802.5, 7: 198.0, 8: 450.5, 9: 0, 10: 416.75, 11: 553.75, 12: 588.25, 13: 448.25, 14: 548.25, 15: 744.0, 16: 295.25, 17: 340.25, 18: 575.25, 19: 594.75, 20: 1000, 21: 1000, 22: 1000, 23: 26.75, 24: 332.0, 25: 640.5, 26: 718.75, 27: 1000, 28: 959.25, 29: 137.75, 30: 15.0, 31: 243.75, 32: 1000, 33: 210.0, 34: 871.5, 35: 585.5, 36: 1000, 37: 1000, 38: 130.0, 39: 614.5, 40: 1000, 41: 1000, 42: 1000, 43: 255.0, 44: 1000, 45: 290.75, 46: 180.25, 47: 530.25, 48: 132.0, 49: 8.25, 50: 520.0}

three= {1: 958.0, 2: 470.0, 3: 0, 4: 0, 5: 827.5, 6: 0, 7: 0, 8: 545.0, 9: 802.5, 10: 425.0, 11: 398.0, 12: 426.75, 13: 198.5, 14: 588.75, 15: 548.25, 16: 742.0, 17: 0, 18: 296.0, 19: 254.5, 20: 337.75, 21: 638.75, 22: 570.25, 23: 447.0, 24: 26.75, 25: 0, 26: 715.75, 27: 210.0, 28: 136.75, 29: 16.0, 30: 246.5, 31: 0, 32: 63.5, 33: 0, 34: 0, 35: 588.0, 36: 425.5, 37: 0, 38: 126.25, 39: 611.25, 40: 719.0, 41: 0, 42: 299.25, 43: 179.5, 44: 130.75, 45: 543.0, 46: 10.5, 47: 0, 48: 0, 49: 0, 50: 0}

four={1: 952.5, 2: 470.0, 3: 1000, 4: 1000, 5: 828.0, 6: 1000, 7: 1000, 8: 545.5, 9: 803.25, 10: 425.5, 11: 398.5, 12: 416.5, 13: 197.0, 14: 586.5, 15: 548.0, 16: 740.25, 17: 1000, 18: 294.5, 19: 254.5, 20: 336.5, 21: 639.25, 22: 570.75, 23: 446.0, 24: 24.0, 25: 1000, 26: 713, 27: 210.25, 28: 137.0, 29: 14.25, 30: 243.75, 31: 1000, 32: 63.75, 33: 1000, 34: 1000, 35: 588.75, 36: 374.75, 37: 1000, 38: 126.75, 39: 610.25, 40: 719.25, 41: 1000, 42: 294.5, 43: 176.75, 44: 131.5, 45: 530.5, 46: 9.0, 47: 1000, 48: 1000, 49: 1000, 50: 1000}
#%%
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import random
import sys

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
atominfo = bfp.parsebondfile(bondfile,bo=True)
#%%--
neighbours = atominfo['neighbours']
atomtypes  = atominfo['atypes']
bondorders = atominfo['bondorders']

steps = list(neighbours.keys())

firststepneigh = neighbours[steps[0]]

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
#%%
plot_keys = ['C-O-upper','C-O-lower','O-H','C-C']

name    = {'O-H':'Bond-1','C-C':'Bond-2','C-O-upper':'Bond-3','C-O-lower':'Bond-4'}
redline = {'O-H':one,'C-C':two,'C-O-lower':three,'C-O-upper':four}

final_temp = 1200
skipts = (final_temp-300)/4
print('skip',skipts)
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Heatmaps'

for key in plot_keys:
    title = '{}-{}\n\nMolecule D'.format(key,name[key])
    timestep = 0.25
    bond = bonds[key]
    
    # getting one, two, three and four
    i = 1
    bo = {}
    for atom1,atom2 in bond:
        bo[i] = min([bondorders[x][atom1].get(atom2,100) for x in steps])
        i+=1
    print(bo)
    sys.exit()
    bfp.bondorder_evolution(bondorders, bond, ts=timestep,savedir=None,title=title,skipts=skipts,figsize=(6,10))
    for i,red in enumerate(bo):
        vline = bo[red]
        plt.plot([vline*4]*10,np.linspace(i, i+1,10),color='red')
        # if vline<50:
        #     plt.text((vline+10)*4,i+0.8,'{}'.format(vline),fontsize=8,fontweight="bold",color='white')
        # else:
        #     plt.text((vline-60)*4,i+0.8,'{}'.format(vline),fontsize=8,fontweight="bold",color='white')
    sd = savedir+'\\bo_evolution_'+str(random.randint(100000, 999999))
    plt.xticks(rotation=0,fontsize=12)
    plt.savefig(sd, dpi=400, bbox_inches='tight')