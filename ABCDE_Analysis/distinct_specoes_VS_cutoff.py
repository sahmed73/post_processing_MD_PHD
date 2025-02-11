# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Mar 27 19:01:18 2024
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import magnolia.speciesfile_parser as sfp
import seaborn as sns

# Set default font sizes
plt.rc('font', size=22) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=22) # Axes label size
plt.rc('xtick', labelsize=20) # X-axis tick label size
plt.rc('ytick', labelsize=20) # Y-axis tick label size
plt.rc('legend', fontsize=15) # Legend fontsize

plt.rcParams['font.size'] = 12

               
labels=['Molecule A','Molecule B']
colors=['tab:red','tab:blue']
fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True, sharex=True)

dirrs = [r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1000K",
        r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1250"]

for i, dirr in enumerate(dirrs):
    n_distinct_species = {}
    sim_dir = ['Sim-1','Sim-2','Sim-3']
    for sim in sim_dir:
        n_distinct_species[sim] = {}
        for cutoff in np.array(range(30,76,5))/100:
            directory = dirr+'\\'+sim
            speciesfile = directory+f'\\species_{cutoff}.out'
            if cutoff==0.30:
                speciesfile = directory+'\\species.out'
                
            species_count = sfp.get_species_count(speciesfile)
            n_distinct_species[sim][cutoff] = species_count.columns.size
            
    df = pd.DataFrame(n_distinct_species)
    
    ax[i].errorbar(df.index, df.mean(axis=1), yerr=df.std(axis=1),
                fmt='-s', capsize=5, ecolor='k', label=labels[i],
                color=colors[i])
    ax[i].tick_params(axis='both', which='major', top=True, bottom=True,
                   left=True, right=True)
    ax[i].set_xticks(np.arange(0.3,0.81,0.1))
        
    # ax[i].legend()
    
    ax[i].set_xlabel('Bond order cutoff')
    ax[i].grid(False)
ax[0].set_ylabel('# distinct species')
ax[0].set_yticks(range(0,450,100))
ax[0].set_ylim(-15,430)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+'\\species_speciesfile{}'.format(random.randint(0,10000000000)),dpi=400,bbox_inches='tight')