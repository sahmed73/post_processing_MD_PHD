# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 13:51:11 2023

@author: Shihab
"""

import magnolia.bondfile_parser as bfp
import magnolia.log_parser_FHB as lfp
import magnolia.needless_essential as ne
import pandas as pd
import matplotlib.pyplot as plt
import os
import random
import time
import numpy as np
import sys
import pickle


start_time = time.time()

timestep = 0.25
print('Assigned Timestep: ',timestep)

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO2\50_AO2_300_O2\Production\Sim-1'

filename    = 'bonds.reaxc'
bondfile    = directory+'\\'+filename

# check if the bond pickle file exists
picklefile     = directory+ '\\' + 'pickle_bond.pkl'
if os.path.exists(picklefile):
    print('Bondfile is loading from Pickle...')
    load_timestep, output_get_neighbours = pickle.load(open(picklefile, 'rb'))
    if load_timestep!=timestep:
        sys.exit('Error: Assigned timestep does not match with loaded the timestep\nAssigned Timestep: {} and Loaded Timestep: {}'.format(timestep,load_timestep))
else:
    print('Calling get_neighbours function....')
    output_get_neighbours = bfp.get_neighbours(bondfile,bo=True)
    with open(picklefile,'wb') as pf:
        pickle.dump((timestep,output_get_neighbours), pf)

neighbours,atomtypes,bondorders = output_get_neighbours

steps = list(neighbours.keys())

print('bondfile is loaded!',end='\t')
ne.print_runtime(time.time()-start_time)

# check if the SpeciesCountAtEveryTimestep pickle file exists
picklefile = directory + '\\'+ 'pickel_SpeciesCountAtEveryTimestep.pkl'
if os.path.exists(picklefile):
    print('SpeciesCountAtEveryTimestep is loading from Pickle...')
    load_timestep, data = pickle.load(open(picklefile, 'rb'))
    if load_timestep!=timestep:
        sys.exit('Error: Assigned timestep does not match with loaded the timestep\nAssigned Timestep: {} and Loaded Timestep: {}'.format(timestep,load_timestep))
else:
    print('Calling get_SpeciesCountAtEveryTimestep function....')
    data = bfp.get_SpeciesCountAtEveryTimestep(neighbours, atomtypes, 'HCO',step2ps=True,step2psargs=[timestep])
    with open(picklefile,'wb') as pf:
        pickle.dump((timestep,data), pf)

print('get_SpeciesCountAtEveryTimestep is loeded!',end='\t')
ne.print_runtime(time.time()-start_time)
#%%
#thermo = lfp.thermo_dict(directory+'\\log.lammps', 2)
#lfp.tempramp(thermo, timestep,Print=True)

df = pd.DataFrame(data)
df = df.fillna(0)

AO_1 = 'H12C10O3' 
AO_2 = 'H34C26O4'
seek_species =  ['H2CO','H4CO','H2O','CO2','O2']
label = bfp.make_molecular_formula_latex(seek_species,sort=True) 
color = ['red','green','blue','black','orange']

for i in range(len(seek_species)):
    if seek_species[i] not in df.index:
        df.loc[seek_species[i],:]=0
        

fig,ax = plt.subplots()

for i in range(len(seek_species)):
    species = df.loc[seek_species[i]]
    ax.plot(species,
            color=color[i],
            label=label[i])
ax.set_xlabel("Time (ps)")
ax.set_ylabel("Number of molecules")
#ax.set_ylim(top=205)
#ax.yaxis.set_ticks(np.arange(0, 201, step=25))
#ax.plot([],[], color='orange',label='$O_2$')
ax.legend(loc='center left')

# ax2=ax.twinx()
# ax2.plot(df.loc['O2',],color='orange',label='$O_2$')
# ax2.set_ylabel("Number of ${O_2}$")
#ax2.set_ylim([45,210])
#ax2.yaxis.set_ticks(np.arange(50, 210, step=50))


ax.set_title(directory[directory.find('LAMMPS'):]+'\n\n')
savedir = 'python_outputs\\figures\\species_plot'+str(random.randint(10000, 99999999))
fig.savefig(savedir, dpi=300, bbox_inches='tight')
plt.show()

print('--------------------------')
ne.print_runtime(time.time()-start_time)
print('--------------------------')