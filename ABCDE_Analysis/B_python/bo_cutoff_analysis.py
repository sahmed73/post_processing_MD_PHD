# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec 29 03:52:55 2023
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import magnolia.speciesfile_parser as sfp
               
data = {}
all_data = {}
total_species = {}

for cutoff in np.array(range(30,80,5))/100:
    total_species[cutoff]=0
    
for cutoff in np.array(range(30,80,5))/100:
    dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1000K"
    sim_dir = ['\\Sim-1','\\Sim-2','\\Sim-3']
    flag = True
    for sim in sim_dir:
        directory = dirr+sim
        speciesfile = directory+f'\\species_{cutoff}.out'
        if cutoff==0.30:
            speciesfile = directory+'\\species.out'
        current = sfp.get_species_count(speciesfile).T
        if flag:
            summed_species = current.copy()
            flag = False
        else:
            summed_species = summed_species.add(current,fill_value=0)
    exclude= []#['O2', 'H26C34O4']
    nsim   = len(sim_dir) # number of simulation
    species = (summed_species/nsim).drop(exclude,axis=1)
    plt.style.use("classic")
    plt.figure(figsize=(6.4, 4.8))
    total_molecule = species.count(axis=1)
    total_species[cutoff]+= len(species.columns)
#%%   
plt.plot(total_species.keys(),total_species.values(),marker='s',c='r')
plt.xlim(0.29)
plt.xlabel('Bond order cutoff')
plt.ylabel('Total number of different species')
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+'\\species_speciesfile{}'.format(random.randint(0,10000000000)),dpi=400,bbox_inches='tight')