# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 01:16:37 2023

@author: arup2
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
import re


start_time = time.time()

timestep = 0.25
print('Assigned Timestep: ',timestep)

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Squalane\Production\1350\Sim-1'
filename    = 'bonds.reaxc'
bondfile    = directory+'\\'+filename

neighbours,atomtypes = bfp.get_neighbours(bondfile)
steps = list(neighbours.keys())

data = bfp.get_SpeciesCountAtEveryTimestep(neighbours, atomtypes, 'HCO',step2ps=True,step2psargs=[timestep])
#%%
#thermo = lfp.thermo_dict(directory+'\\log.lammps', 2)
#lfp.tempramp(thermo, timestep,Print=True)

df = pd.DataFrame(data)
df = df.fillna(0)
lignins = {'A':'H12C10O3',
           'B':'H34C26O4',
           'C':'H44C29O2',
           'D':'H30C19O3',
           'E':'H14C11O4',
           }

#seek_species =  ['H2CO','H4CO','H2O','CO2','H4CO2']
pattern = r'\b\d{3,4}\b'
matches = re.findall(pattern, directory)
print(matches)
final_temp = int(matches[0])
skipts = (final_temp-300)/4
seek_species = ['O2','H62C30']
label = bfp.make_molecular_formula_latex(seek_species,sort=True) 
color = ['red','#66CCFF','blue','black','maroon','orange']
marker = ['s','^','v','o','*','>']

for i in range(len(seek_species)):
    if seek_species[i] not in df.index:
        df.loc[seek_species[i],:]=0
        
df = df.loc[:,skipts:]
df.columns = df.columns-skipts       

fig, ax = plt.subplots()

for i in range(len(seek_species)):
    species = df.loc[seek_species[i],:]
    ps = np.array(species.index)
    count = np.array(species.values)
    if label[i]!='$O_{2}$':
        lab = 'Squalane'
    else:
        lab = label[i]
    ax.plot(ps, count,
            color=color[i],
            label=lab)
ax.set_xlabel("Time (ps)")
ax.set_ylabel("Number of molecules")
ax.legend()
ax.set_xlim(0,1000+10)
ax.set_ylim(bottom=0)
#ax.yaxis.set_ticks(np.arange(0, 201, step=25))

ax.set_title(directory[directory.find('LAMMPS'):]+'\n\n'+'Molecule: D')
savedir = '..\\python_outputs\\figures\\species_plot'+str(random.randint(10000, 99999999))
fig.savefig(savedir, dpi=300, bbox_inches='tight')
plt.show()

print('--------------------------')
ne.print_runtime(time.time()-start_time)
print('--------------------------')