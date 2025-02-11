# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 23 05:15:03 2023
"""

import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
import random
import pandas as pd
import numpy as np

base_oil = 'PAO4'
    
commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600".format(base_oil,base_oil)

for bo_cutoff in np.array(range(35,76,5))/100:
    flag = True
    simulation_dirs = ['\\Sim-1','\\Sim-2','\\Sim-3']
    for sim in simulation_dirs:#,'\\Sim-2','\\Sim-3']:
        print(sim)
        directory = commmon+sim
        bondfile = directory+'\\bonds.reaxc'
        current = bfp.get_species_count(bondfile,['H','C','O'],cutoff=bo_cutoff).T
        if flag:
            summed_species = current.copy()
            upper_bound_   = current.copy()
            lower_bound_   = current.copy()
            flag = False
        else:
            summed_species = summed_species.add(current,fill_value=0)
            ## geting upper bound
            maxx = lambda s1, s2: s1.where(s1 > s2, s2)
            upper_bound_ = upper_bound_.combine(current, maxx)
            ## getting lower bound
            minn = lambda s1, s2: s1.where(s1 < s2, s2)
            lower_bound_ = lower_bound_.combine(current, minn)
    
    #%%
    number_of_main = 25  
    cutoff = 1
    skipts = (1600-300)/4
    exclude= ['O2', 'H62C30']
    nsim   = len(simulation_dirs) # number of simulation
    
    species     = summed_species/nsim
    upper_bound = upper_bound_.copy()
    lower_bound = lower_bound_.copy()
    
    last_count = species.iloc[-1,:].sort_values(ascending=False)
    species    = species.loc[:,last_count.index]
    species = species.loc[:,(species>=cutoff).any()]
    species = species.drop(exclude,axis=1)
    species.index = species.index*0.25/1000
    species = species.loc[skipts:,:]
    species.index = species.index-skipts
    
    ## adjusting the index of upper_bound
    upper_bound.index = upper_bound_.index*0.25/1000
    upper_bound = upper_bound.loc[skipts:,:]
    upper_bound.index = upper_bound.index-skipts
    ## adjusting the index of lower_bound
    lower_bound.index = lower_bound_.index*0.25/1000
    lower_bound = lower_bound.loc[skipts:,:]
    lower_bound.index = lower_bound.index-skipts
    
    
    time = species.index
    
    colors = iter(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    plt.style.use("classic")
    plt.figure(figsize=(6.4, 4.8))
    
    for column in species.columns:
        label = bfp.make_molecular_formula_latex(column,sort=True)
        # color = next(colors)
        plt.plot(time,species[column],label=label,linewidth=0.8)
        plt.fill_between(time, lower_bound[column], upper_bound[column],
                         alpha=0.1)
        
    plt.legend(fontsize=7,loc='upper left')
    plt.ylim(0,60)
    plt.xlabel('Time (ps)',fontsize=15)
    plt.ylabel('Number of molecules',fontsize=15)
    title = f'This plot is generated from bond file with bond order cutoff {bo_cutoff}'
    plt.title(directory[directory.find('Base'):]+'\n'+title+'\n',fontsize=8)
    plt.grid('on')
    plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+'\\species_bondfile_{}'.format(random.randint(0,10000000000)),dpi=400,bbox_inches='tight')
    plt.show()