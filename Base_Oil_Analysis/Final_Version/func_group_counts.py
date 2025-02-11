# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Oct 31 15:30:38 2023
"""

import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
import random
import numpy as np

### Directory ###
base_oil = 'Squalane'
common = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600".format(base_oil,base_oil)

fg_counts = []

for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    directory = common+sim
    filename  = "\\bonds.reaxc"
    bondfilepath = directory+filename
    
    ### Parsing Bondfile ###
    bonddata = bfp.parsebondfile(bondfilepath)
    
    neighbours = bonddata['neighbours']
    atypes     = bonddata['atypes']
    asyms      = ['H','C','O']
        
    
    ##searching functional groups
    seek = ['OH','COOH','Keto', 'Aldy']
    name = ['Hydroxyl','Carboxylic acid', 'Ketone', 'Aldehyde']
    
    count = bfp.count_functinalGroup(neighbours,atypes,seek)
    fg_counts.append(count)
#%%
fg_counts_array = np.array(fg_counts)
means = fg_counts_array.mean(axis=0)
stds   = fg_counts_array.std(axis=0)



plt.style.use("classic")
plt.figure(figsize=[6,4.5])
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\functional_group\\fg_barplot_{}'.format(random.randint(0,999999999))
title = directory[directory.find('Base'):]+'\n'
plt.title(title,fontsize=7)
plt.bar(name,means,yerr=stds, align='center',alpha=0.9, ecolor='black',
        capsize=10,color='maroon')
plt.ylabel('Number of molecules')
plt.ylim(0,160)
plt.savefig(savedir, dpi=400, bbox_inches='tight')