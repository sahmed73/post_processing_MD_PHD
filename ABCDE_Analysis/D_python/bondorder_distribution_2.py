# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Apr 23 10:59:15 2024
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import magnolia.plot_template as mplt
from matplotlib.ticker import AutoMinorLocator, LogLocator

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1000K\Sim-1'
filename  = '\\bonds.reaxc'
bondfile_D  = dirr+filename
bonddata_D  = bfp.parsebondfile(bondfile_D,bo=True)
#%%
neighbors = bonddata_D['neighbours']
atypes    = bonddata_D['atypes']
bondorders= bonddata_D['bondorders']

bolist = []
for step in bondorders:
    for u in bondorders[step]:
        for v in bondorders[step][u]:
            bo = bondorders[step][u][v]
            bolist.append(bo)
boarrD = np.array(bolist)
#%%
dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1250\Sim-1'
filename  = '\\bonds.reaxc'
bondfile_B  = dirr+filename
bonddata_B  = bfp.parsebondfile(bondfile_B,bo=True)
#%%
neighbors = bonddata_B['neighbours']
atypes    = bonddata_B['atypes']
bondorders= bonddata_B['bondorders']

bolist = []
for step in bondorders:
    for u in bondorders[step]:
        for v in bondorders[step][u]:
            bo = bondorders[step][u][v]
            bolist.append(bo)
boarrB = np.array(bolist)
#%%
savedir = (r'C:\Users\arup2\OneDrive - University of California Merced'
           r'\Desktop\LAMMPS\Post Processing\lammps_processing'
           r'\python_outputs\figures')

# Set default font sizes
plt.rcdefaults()
plt.rc('font', size=22) # Default text sizes
plt.rc('axes', titlesize=22) # Axes title size
plt.rc('axes', labelsize=22) # Axes label size
plt.rc('xtick', labelsize=45) # X-axis tick label size
plt.rc('ytick', labelsize=45) # Y-axis tick label size
plt.rc('legend', fontsize=15) # Legend fontsize

labels=['Molecule A','Molecule B']
colors=['tab:red','tab:blue']
fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=0, sharex=False)
data = [boarrD,boarrB]
fig.subplots_adjust(wspace=0.7)
for i in range(2):
    axs[i].hist(data[i], bins=80, color=colors[i], density=True, 
             edgecolor='none', label=labels[i])
    # axs[i].set_yscale('log')
    axs[i].tick_params(axis='both', which='major', top=True, bottom=True,
                   left=True, right=True)
    # axs[i].set_xticks(np.arange(0.0,3.1,0.5))
    # axs[i].set_yticks(np.arange(0,13,2))
    axs[i].grid()
    # axs[i].set_title(labels[i]+'\n')
    # axs[i].set_xlabel('Bond order')
    axs[i].set_xlim(0.3,0.5)
    axs[i].grid(False)
    axs[i].tick_params(axis='x', which='major', pad=15)
    axs[i].set_ylim(0,0.004)
    # axs[i].set_ylim(0.0,1e7)
    # axs[i].yaxis.set_minor_locator(
    #     plt.LogLocator(base=10.0,subs=np.arange(2, 10) * 0.1, numticks=10))
    # axs[i].xaxis.set_minor_locator(AutoMinorLocator())
    

# fig.text(0.5, -0.05, 'Bond order', ha='center')
# axs[0].set_ylabel('Frequency')
# axs[0].set_yticks(range(2,13,2))
mplt.saveplot('figures', 'bo_dist')

