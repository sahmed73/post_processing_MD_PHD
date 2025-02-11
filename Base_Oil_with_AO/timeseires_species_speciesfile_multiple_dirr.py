# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Mar 24 04:34:16 2024
"""

import magnolia.speciesfile_parser as sfp
import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
import random
import pandas as pd
import numpy as np
import magnolia.plot_template as mplt

# Set default font sizes
plt.rc('font', size=15) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=18) # Axes label size
plt.rc('xtick', labelsize=15) # X-axis tick label size
plt.rc('ytick', labelsize=15) # Y-axis tick label size
plt.rc('legend', fontsize=12) # Legend fontsize

base_oil = 'PAO4'
molecules = ['H30C29O3', 'H34C26O4']
molecule = molecules[1]

temps = ['1250K','1150K','1050K']
    
locations = [r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1250",
            r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1150",
            r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1050"]

fig, ax = plt.subplots()

for i, location in enumerate(locations):

    sim_dir = ['Sim-1', 'Sim-2','Sim-3']
    species = {}
    dfs = pd.DataFrame()
    for sim in sim_dir:
        ####################################
        ############-ALERT--##########################
        if i==2 and sim!='Sim-1': continue
        ########################################
        
        
        speciesfilepath = location+"\\"+sim+'\\species.out'
        species[sim] = sfp.get_species_count(speciesfilepath)
        df = species[sim].T[molecule]
    
        timestep        = 0.25
        ramp_rate       = 4
        initial_temp    = 300
        skipts          = 0#(900-initial_temp)/ramp_rate
        
        df.index = df.index*timestep/1000-skipts
        df = df.loc[0:]
        
        dfs[sim] = df.copy()
    
    mplt.mean_range_plot(dfs, label=temps[i], ax=ax, shade=False)
    
saveplt = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures'+f'\\timeseries_species_{random.randint(0,100000000)}'

ax.legend()
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Number of molecule B')
fig.savefig(saveplt,dpi=300,bbox_inches='tight')