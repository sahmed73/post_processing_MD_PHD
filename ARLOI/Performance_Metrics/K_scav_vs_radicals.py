# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri May 24 02:08:26 2024
"""
import sys
import re
import magnolia.speciesfile_parser as sfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('default')
plt.rcParams['font.size'] = 18

molecule='B'
dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\20PAOtcr_15{}\Production\TempRamp'.format(molecule)
filename = '\\species.out'
speciesfile = dirr+filename
species = sfp.get_species_count(speciesfile, timestep=0.25, Nfreq = 1)

radicals = species['H61C30O']
time = species['Time']/1000 # ns
# print(species.columns)
# pattern = re.compile(r'O(?!\d)')
# s = [x for x in species.columns if pattern.search(x)]
# if 'H62C30O' in s: s.remove("H62C30O")
# print(s)
# print(species[s].sum(axis=1))
# sys.exit(0)

## data stripping
N_values = [8]
m_values = []
fig, ax = plt.subplots(dpi=350)
for N_scav in N_values:
    N = radicals.max()-N_scav
    print(N,N_scav)
    idx=None if N is None else (radicals[::-1]==N).idxmax()+1
    radicals=radicals.loc[:idx]
    time=time.loc[:idx]
    log_radicals = np.log(radicals)
    
    ## fitting
    m,c  = np.polyfit(time, log_radicals, 1)
    m_values.append(-m)
    
    # Convert to scientific notation
    base, exp = f"{-m:.2e}".split('e')
    k_scav = f"{base}x10^{{{int(exp)}}}"
    
    
    # plotting
    ax.scatter(time, log_radicals, s=4, alpha=0.8) #, label='data')
    label=f"$N_{{scav}}={N_scav}, k_{{scav}}={-m:.4} \ ns^{{-1}}$"
    ax.plot(time,m*time+c, label=label, color='k')
ax.set_xlabel("Time (ns)")
ax.set_ylabel("ln(N$_{radical}$)")
ax.set_title(f"Antioxidant {molecule}", fontsize=16)
ax.legend(fontsize=12)
ax.set_ylim(1.7, 2.5)

# fig, ax = plt.subplots(dpi=350)
# ax.plot(N_values,m_values,'o-')
# ax.set_xlabel("Number of scavanged radical")
# ax.set_ylabel("$k_{scav} \ (ns^{-1})$")
# ax.set_title(f"Antioxidant {molecule}", fontsize=16)
# # ax.legend(fontsize=10)
# ax.set_ylim(0.8, 3.3)

# ax.text(0.45,0.9,f"$k_{{scav}}={-m:.4} \ ns^{{-1}}$",transform=ax.transAxes, fontsize=16)