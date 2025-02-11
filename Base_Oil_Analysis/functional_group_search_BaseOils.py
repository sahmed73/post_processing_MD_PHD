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
base_oil = 'PAO4'
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(base_oil,base_oil)
filename  = "\\bonds.reaxc"
bondfilepath = directory+filename

### Parsing Bondfile ###
bonddata = bfp.parsebondfile(bondfilepath)
#%% Counting functional group
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
atomsymbols= 'HCO'

##searching functional groups
seek = ['pp','OH','COOH','Keto', 'Aldy']
name = ['Hydroxyl','Carboxylic acid', 'Ketone', 'Aldehyde']

group = {'OH': [3,[1,2]],
         'COOH': [2,[1,3,3]],
         'Keto': [2,[2,2,3]],
         'Aldy': [2,[1,2,3]]}

count = bfp.count_functinalGroup(neighbours, atypes, seek)
#%%
plt.style.use("classic")
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\functional_group\\fg_barplot_{}'.format(random.randint(0,999999999))
title = directory[directory.find('Base'):]+'\n'
plt.title(title,fontsize=7)
plt.bar(name,count,color='maroon',alpha=0.9)
plt.ylabel('Number of molecules')
plt.ylim(0,160)
plt.savefig(savedir, dpi=400, bbox_inches='tight')