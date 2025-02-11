# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon May  6 03:04:50 2024
"""

import magnolia.speciesfile_parser as sfp
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt
import scienceplots
import magnolia.plot_template as mplt
plt.style.use('default')
# plt.style.use('science')
plt.rcParams['font.size'] = 18

directories = [r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3A_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K',
               r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3B_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K',
               r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3BHT_Using_PACKMOL\Production\StartedFromEnergyEquilibration\600K']

end = 2000
seek      = ['H30C19O3', 'H34C26O4', 'H24C15O']
seek_name = ['A', 'B', "BHT"] 

fig, axs = plt.subplots(1,len(directories),dpi=300, sharey=True)
for i, (ax, dirr) in enumerate(zip(axs,directories)):   
    filename = '\\species.out'
    speciesfile = dirr+filename
    species = sfp.get_species_count(speciesfile, timestep=0.25)

    time = species['Time'].iloc[:end]
    species = species.iloc[:,:-1]
    ax.plot(time, species['H62C30O'].iloc[:end])
    ax.set_title(seek_name[i])
    
axs[0].set_ylabel('Number of SR')
fig.text(0.4, -0.05, 'Time (ps)')