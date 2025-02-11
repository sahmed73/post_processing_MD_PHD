# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu May 23 06:55:04 2024
"""

import magnolia.speciesfile_parser as sfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('default')
plt.rcParams['font.size'] = 18

rad_name = 'A'

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\research2024\Temp\20_PAO-OH_15_BHT\Production\Ashraf_2016'
filename = '\\species.out'
speciesfile = dirr+filename
species = sfp.get_species_count(speciesfile, timestep=0.25, Nfreq = 1)

radicals = species['H61C30O']
time = species['Time']/1000 # ns

### scaling
# n_radicals = 20
# missing = n_radicals - radicals.max()
# r = radicals/radicals.max()
# radicals+= r*missing


## data stripping
N=4
idx=None if N is None else (radicals[::-1]==N).idxmax()+1
radicals=radicals.iloc[:idx]
time=time.iloc[:idx]
log_radicals = np.log(radicals+0.001)


## fitting
m,c  = np.polyfit(time, log_radicals, 1)

# Convert to scientific notation
base, exp = f"{-m:.2e}".split('e')
k_scav = f"{base}x10^{{{int(exp)}}}"


# plotting
fig, ax = plt.subplots(dpi=350)
ax.scatter(time, log_radicals, s=4, alpha=0.8, label='data')
ax.plot(time,m*time+c,color='k', label='fit')
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Number of radicals")
ax.set_title(f"Molecule {rad_name}", fontsize=16)
ax.set_xlim(-0.025, 0.425)
ax.set_ylim(1-0.05,3.2)
# ax.legend()
ax.text(0.47,0.92,f"$k_{{scav}}={-m:.4} \ ns^{{-1}}$",transform=ax.transAxes,
        fontsize=16)