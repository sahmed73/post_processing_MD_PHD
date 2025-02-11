# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Feb  1 03:31:13 2024
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

bonddata = {}
cutoff   = 0.35
for sim in simdir:
    directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(baseoil,baseoil,sim)
    bondfilepath = directory+'\\bonds.reaxc'
    bonddata[sim] = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%%
sim = 'Sim-1'

neighbors = bonddata[sim]['neighbours']
atypes    = bonddata[sim]['atypes']
mtypes    = bonddata[sim]['mtypes']
bondorders= bonddata[sim]['bondorders']
atomsymbols = 'HCO'

atom_type_pairs = [[1,2],[1,3],[2,2],[2,3],[3,3]]
atom_symb_pairs = ['C-H','O-H','C-C','C-O','O-O']
pairs = zip(atom_type_pairs,atom_symb_pairs)

for check, symb in pairs:
    bond_order = []
    print(check,symb)
    for step in bondorders:
        for u in bondorders[step]:
            for v in bondorders[step][u]:
                pair = sorted([atypes[u],atypes[v]])
                bo = bondorders[step][u][v]
                if check==pair and bo not in bond_order:
                    bond_order.append(bo)

    bins = 20
    label_fontsize = 15
    savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\bond_order_hist_{}_{}_bins_{}'.format(symb,bins,random.randint(0,100000000))
    fig, ax = plt.subplots()
    counts, edges, patches = ax.hist(bond_order,bins=bins,
                                     color='tab:red')
    
    # Adding text on top of each bar
    total_count = sum(counts)

# Adding percentage text on top of each bar
    for i in range(len(counts)):
        percentage = 20 * counts[i] / total_count
        # plt.text(edges[i] + (edges[1] - edges[0]) / 2, counts[i], f'{percentage:.1f}%', ha='center', va='bottom', fontsize=7, rotation=45)

        
        
        
    hist_min = edges[0]
    hist_max = edges[-1]
    print(counts)
    print(edges)
    print('--------------------------------------')
    
    text_properties = {
        'facecolor': 'white',  # Box background color
        'edgecolor': 'black',  # Box edge color
        'boxstyle': 'round,pad=0.5'  # Box style and padding
    }
    
    ax.text(0.95, 0.95, f'{symb} Bonds\nMin BO={hist_min}\nMax BO={hist_max}\n{bins} bins', horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, bbox=text_properties)
    ax.set_xlabel('Bond ordes',fontsize=label_fontsize)
    ax.set_ylabel('Counts', fontsize=label_fontsize)
    xticks = np.array(range(30,276,30))/100
    ax.set_xlim(xticks[0]-0.1,xticks[-1]+0.1)
    ax.set_xticks(xticks)
    fig.savefig(savedir, dpi=500, bbox_inches='tight')
    plt.show()