# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon May  6 12:52:09 2024
"""

import magnolia.speciesfile_parser as sfp
import magnolia.bondfile_parser as bfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt
plt.style.use('default')
plt.rcParams['font.size'] = 16

directories = [r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3A_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K',
               r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3B_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K',
               r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3BHT_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K']

fig, axs = plt.subplots(1,3, sharey=True, figsize=[5,3], dpi=300)
fig.subplots_adjust(wspace=0.5)
seek      = ['H30C19O3', 'H34C26O4', 'H24C15O']
seek_name = ['A', 'B', "BHT"] 

for k, dirr in enumerate(directories):
    filename = '\\species.out'
    speciesfile = dirr+filename
    species = sfp.get_species_count(speciesfile, timestep=0.25)
    time = species['Time']
    
    AO  = species[seek[k]] 
    RAD = species['H62C30O']
    if k==1:
        RAD = species['H62C30O']
    
    axs[k].plot(time, AO)
    axs[k].plot(time,RAD)
    axs[k].set_ylim(0,3)
    axs[k].set_title(seek_name[k])


fig.text(0.5, -0.06, 'Time (ns)', ha='center')
axs[0].set_ylabel('Number of molecule')
mplt.saveplot('species','species')