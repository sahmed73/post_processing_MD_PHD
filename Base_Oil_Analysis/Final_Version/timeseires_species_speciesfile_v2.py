# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Jan  6 01:54:31 2024
"""

import pandas as pd
import matplotlib.pyplot as plt
import random
import magnolia.speciesfile_parser as sfp
import magnolia.bondfile_parser as bfp
import magnolia.plot_template as mplt
               
                

base_oil = 'PAO4'
cutoff   = 0.3
temp     = 1600
    
commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\{}".format(base_oil,base_oil,temp)

sim_dir = ['\\Sim-1','\\Sim-2','\\Sim-3']
flag = True
all_species = {}
for sim in sim_dir:
    directory = commmon+sim
    speciesfile = directory+f'\\species_{cutoff}.out'
    if cutoff==0.30:
        speciesfile = directory+'\\species.out'
        
    all_species[sim] = sfp.get_species_count(speciesfile).T





number_of_main = 25  
species_cutoff = 4
skipts = (temp-300)/4
exclude= ['O2', 'H62C30']
columns = ['H2O','O2','H62C30']

fig, ax = plt.subplots()
for column in columns:
    mplt.mean_range_plot(all_species,column=column,ax=ax)
    
plt.legend(fontsize=15,loc='upper left')
# plt.ylim(0,60)
plt.xlabel('Time (ps)',fontsize=15)
plt.ylabel('Number of molecules',fontsize=15)
title = f'This plot is generated from species file (cutoff: {cutoff})'
plt.title(directory[directory.find('Base'):]+'\n'+title+'\n',fontsize=8)
plt.grid('on')
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+'\\species_speciesfile{}'.format(random.randint(0,10000000000)),dpi=400,bbox_inches='tight')