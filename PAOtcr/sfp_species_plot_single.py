# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun May  5 22:50:15 2024
"""

import magnolia.speciesfile_parser as sfp
import magnolia.bondfile_parser as bfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt
import scienceplots
plt.style.use('default')
plt.style.use('science')
plt.rcParams['font.size'] = 12

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3BHT_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K'
filename = '\\species.out'
speciesfile = dirr+filename
species = sfp.get_species_count(speciesfile, timestep=0.25)
print(species)

fig, ax = plt.subplots(dpi=300)
time = species['Time']
species = species.iloc[:,:-1]
ax.plot(time,species['H62C30O'])

# for spec in species:
#     print(spec)
#     ax.bar(time,species[spec])