# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec  8 13:29:43 2023
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import magnolia.speciesfile_parser as sfp
               
base_oil = 'PAO4'
data = {}
all_data = {}
total_species = {}
for cutoff in np.array(range(30,76,5))/100:
    commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600".format(base_oil,base_oil)
    sim_dir = ['\\Sim-1','\\Sim-2','\\Sim-3']
    flag = True
    for sim in sim_dir:
        directory = commmon+sim
        speciesfile = directory+f'\\species_{cutoff}.out'
        if cutoff==0.30:
            speciesfile = directory+'\\species.out'
        current = sfp.get_species_count(speciesfile).T
        if flag:
            summed_species = current.copy()
            flag = False
        else:
            summed_species = summed_species.add(current,fill_value=0)
    exclude= ['O2', 'H62C30']
    nsim   = len(sim_dir) # number of simulation
    species = (summed_species/nsim).drop(exclude,axis=1)
    plt.style.use("classic")
    plt.figure(figsize=(6.4, 4.8))
    total_molecule = species.count(axis=1)
    total_species[cutoff] = len(species.columns)
    last_value = total_molecule.iloc[-10:].mean()
    data[cutoff]=last_value
    all_data[cutoff]=species
#%%   
plt.plot(total_species.keys(),total_species.values(),marker='s',c='r')
plt.xlim(0.29)
plt.xlabel('Bond order cutoff')
plt.ylabel('Total number of different species')
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+'\\species_speciesfile{}'.format(random.randint(0,10000000000)),dpi=400,bbox_inches='tight')