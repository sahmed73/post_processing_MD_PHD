# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 23 05:15:03 2023
"""
import pandas as pd
import matplotlib.pyplot as plt
import random
import magnolia.speciesfile_parser as sfp
import magnolia.bondfile_parser as bfp
               
                

base_oil = 'Squalane'
cutoff   = 0.45
temp     = 1600
    
commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\{}".format(base_oil,base_oil,temp)

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





number_of_main = 25  
species_cutoff = 4
skipts = (temp-300)/4
exclude= ['O2', 'H62C30']
nsim   = len(sim_dir) # number of simulation
tailps = 10
nspecies = 5

species     = summed_species/nsim
upper_bound = upper_bound_.copy()
lower_bound = lower_bound_.copy()

species.index = species.index*0.25/1000
species = species.loc[skipts:,:]
species.index = species.index-skipts

last_count = species.loc[1000-tailps:,:].mean().sort_values(ascending=False)
species    = species.loc[:,last_count.index]
# species = species.loc[:,(species>=species_cutoff).any()]
species = species.drop(exclude,axis=1)
species = species.T.head(nspecies).T


## adjusting the index of upper_bound
upper_bound.index = upper_bound_.index*0.25/1000
upper_bound = upper_bound.loc[skipts:,:]
upper_bound.index = upper_bound.index-skipts
## adjusting the index of lower_bound
lower_bound.index = lower_bound_.index*0.25/1000
lower_bound = lower_bound.loc[skipts:,:]
lower_bound.index = lower_bound.index-skipts


## plot parameter
label_fontsize  = 18
tick_fontsize   = 14
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

plt.style.use("default")
fig, ax = plt.subplots(figsize=[8,6])
time = species.index

for i, column in enumerate(species.columns):
    label = bfp.make_molecular_formula_latex(column,sort=True)
    # color = next(colors)
    ax.plot(time,species[column],label=label,linewidth=0.8,
             color=colors[i%len(colors)])
    ax.fill_between(time, lower_bound[column], upper_bound[column],
                     alpha=0.1, color=colors[i%len(colors)])
    
ax.legend(fontsize=tick_fontsize,loc='upper left')
ax.tick_params(axis='both', labelsize=tick_fontsize)
# plt.ylim(0,60)
ax.set_xlim(-1,1001)
ax.set_xlabel('Time (ps)',fontsize=label_fontsize)
ax.set_ylabel('Number of species',fontsize=label_fontsize)
title = f'This plot is generated from species file (cutoff: {cutoff})'
ax.set_title(directory[directory.find('Base'):]+'\n'+title+'\n',fontsize=8)
ax.grid('on')
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+'\\species_speciesfile{}'.format(random.randint(0,10000000000)),dpi=500,bbox_inches='tight')