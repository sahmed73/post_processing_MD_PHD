# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon May  6 01:49:02 2024
"""

import scienceplots
import magnolia.speciesfile_parser as sfp
import matplotlib.pyplot as plt

plt.style.use('default')
# plt.style.use('science')
# plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 18

directories = [r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3A_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K',
               r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3B_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K',
               r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3BHT_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K']


fig, ax = plt.subplots(dpi=500)
seek      = ['H30C19O3', 'H34C26O4', 'H24C15O']
seek_name = ['A', 'B', "BHT"] 
colors    = ['tab:blue','tab:red','tab:green']
times = []

for k, dirr in enumerate(directories):
    filename = '\\species.out'
    speciesfile = dirr+filename
    species = sfp.get_species_count(speciesfile, timestep=0.25)
    time = species['Time']
    
    counts = species[seek[k]]
    for i in range(counts.size):
        if counts.iloc[i]<3:
            print(time.iloc[i])
            times.append(time.iloc[i])
            break

ax.bar(seek_name, times)
ax.set_ylabel('First HAT reaction time (ps)')