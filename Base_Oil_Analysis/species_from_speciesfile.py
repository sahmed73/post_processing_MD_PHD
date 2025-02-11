# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:24:40 2023

@author: arup2
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import magnolia.bondfile_parser as bfp
import random

def get_species(speciesfile):
    with open(speciesfile,'r') as sf:
        species = {}
        for line in sf:
            if "Timestep" in line:
                headers       = line.strip().split()[1:] # not taking '#'
                species_name  = headers[3:] # skip Timestep, No_moles, No_Specs
            else:
                values        = [int(x) for x in line.strip().split()]
                timestep      = values[0]
                species_count = values[3:]
                species[timestep]={}
                for key,value in zip(species_name,species_count):
                    if key in species:
                        species[timestep][key]+=value
                    else:
                        species[timestep][key] = value
    
    df = pd.DataFrame(species).fillna(0).T
    df.index.name = 'Timestep'
    return df
                
                

base_oil = 'Squalane'
    
commmon = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600".format(base_oil,base_oil)

flag = True
for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    directory = commmon+sim
    speciesfile = directory+'\\species.out'
    current = get_species(speciesfile)
    if flag:
        summed_species = current.copy()
        upper_bound    = current.copy()
        lower_bound    = current.copy()
        flag = False
    else:
        summed_species = summed_species.add(current,fill_value=0)
        ## geting upper bound
        maxx = lambda s1, s2: s1.where(s1 > s2, s2)
        upper_bound = upper_bound.combine(current, maxx)
        ## getting lower bound
        minn = lambda s1, s2: s1.where(s1 < s2, s2)
        lower_bound = lower_bound.combine(current, minn)





number_of_main = 25  
cutoff = 4
skipts = (1600-300)/4
exclude= ['O2', 'H62C30']
nsim   = 3 # number of simulation

species = summed_species/nsim

last_count = species.iloc[-1,:].sort_values(ascending=False)
species    = species.loc[:,last_count.index]
species = species.loc[:,(species>=cutoff).any()]
species = species.drop(exclude,axis=1)
species.index = species.index*0.25/1000
species = species.loc[skipts:,:]
species.index = species.index-skipts

## sort upper bound according to the species
upper_bound = upper_bound[species.columns]
upper_bound.index = upper_bound.index*0.25/1000
upper_bound = upper_bound.loc[skipts:,:]
upper_bound.index = upper_bound.index-skipts
# sort lower bound according to the species
lower_bound = lower_bound[species.columns]
lower_bound.index = lower_bound.index*0.25/1000
lower_bound = lower_bound.loc[skipts:,:]
lower_bound.index = lower_bound.index-skipts


time = species.index

colors = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
plt.style.use("classic")
plt.figure(figsize=(6.4, 4.8))

for column in species.columns:
    label = bfp.make_molecular_formula_latex([column],sort=True)
    color = next(colors)
    plt.plot(time,species[column],label=label,linewidth=0.8,c=color)
    plt.fill_between(time, lower_bound[column], upper_bound[column],
                     alpha=0.1,color=color)
    
plt.legend(fontsize=15,loc='upper left')
plt.ylim(0,60)
plt.xlabel('Time (ps)',fontsize=15)
plt.ylabel('Number of molecules',fontsize=15)
plt.title(directory[directory.find('Base'):]+'\n',fontsize=8)
plt.grid('on')
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+'\\species_{}'.format(random.randint(0,10000000000)),dpi=400,bbox_inches='tight')