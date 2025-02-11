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
plt.rcParams['font.size'] = 24

rad_name = ''

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\1_ARLOI_24-08-28--00-09-01\15_PAO-OH_20_A0001\Production\Sim-2_with_H_removed'
filename = '\\species.out'
speciesfile = dirr+filename
species = sfp.get_species_count(speciesfile, timestep=0.25, Nfreq = 1)

radicals = species['H61C30O']
time = species['Time']/1000 # ns

end_number = radicals.iloc[0]*(1-1.0)

time     = time[radicals>end_number]
radicals = radicals[radicals>end_number]

## data stripping
log_radicals = np.log(radicals)


## fitting
m,c  = np.polyfit(time[radicals<20], log_radicals[radicals<20], 1)

# Convert to scientific notation
# base, exp = f"{-m:.2e}".split('e')
# k_scav = f"{base}x10^{{{int(exp)}}}"


# plotting
fig, ax = plt.subplots(dpi=350)
ax.scatter(time, log_radicals, s=4, alpha=0.8, label='data')
ax.plot(time[radicals<20],m*time[radicals<20]+c,color='k', label='fit')
ax.set_xlabel("Time (ns)")
ax.set_ylabel("ln(Number of radicals)")
ax.set_title(f"{rad_name}")
ax.set_xlim(-0.025, 0.425)
ax.set_ylim(-0.1,3.5)
# ax.legend()
ax.text(0.47,0.92,f"$k_{{scav}}={-m:.4} \ ns^{{-1}}$",transform=ax.transAxes, fontsize=16)