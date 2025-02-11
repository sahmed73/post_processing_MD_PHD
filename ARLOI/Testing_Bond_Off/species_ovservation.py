# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Jun  5 09:11:07 2024
"""

import magnolia.bondfile_parser as bfp
import magnolia.speciesfile_parser as sfp
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots(dpi=300)

labels = ['with A', 'with BHT', 'with B']

for i in range(1, 4):
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-{}'.format(i)
    
    file='\\species.out'
    
    species = sfp.get_species_count(dirr+file, timestep=0.25)
    rad     = species.iloc[:,0]*20/25
    time    = species['Time']
    time    = time*425/time.max()
    ax1.plot(time, rad, label=labels[i-1])

ax1.set_xlabel("Time (ps)")
ax1.set_ylabel("Number of radicals")

# Add secondary x-axis on the top for temperature
# secax = ax1.secondary_xaxis('top', functions=(lambda x: 300 + x * 4,
#                                               lambda x: (x-300) / 4))
# secax.set_xlabel("Temperature (K)")
ax1.legend()

plt.show()
