# -*- coding: utf-8 -*-
"""
Created on Fri May 26 07:03:35 2023

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


start_time = time.time()

timestep = 0.25
print('Assigned Timestep: ',timestep)

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\E\E_300_O2\Onset\Sim-1'
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
left = 50

lignins = {'A':'H12C10O3',
           'B':'H34C26O4',
           'C':'H44C29O2',
           'D':'H30C19O3',
           'E':'H14C11O4'}

seek_species =  ['H2CO','H4CO','H2O','CO2','H4CO2']
label = bfp.make_molecular_formula_latex(seek_species,sort=True) 
color = ['red','green','blue','black','maroon','orange']
marker = ['s','^','v','o','*']

for i in range(len(seek_species)):
    if seek_species[i] not in df.index:
        df.loc[seek_species[i],:]=0

print('Sliced')
df = df.loc[:,left:]
df.columns = df.columns - left        

fig,ax = plt.subplots()

for i in range(len(seek_species)):
    species = df.loc[seek_species[i],150:400]
    ax.plot(species,
            color=color[i],
            label=label[i],
            marker=marker[i],
            markevery=100)
ax.set_xlabel("Time (ps)")
ax.set_ylabel("Number of molecules")
#ax.set_xlim(0,1000+10)
#ax.set_ylim([0,350])
#ax.yaxis.set_ticks(np.arange(0, 201, step=25))
ax.plot([],[], color='orange',label='$O_2$',marker='*')
legend = plt.legend(loc="center left", edgecolor="black")
legend.get_frame().set_alpha(None)
legend.get_frame().set_facecolor((0, 0, 1, 0.1))

ax2=ax.twinx()
ax2.plot(df.loc['O2',150:400]/1.25+60,color='orange',label='$O_2$',marker='*',markevery=100)
ax2.set_ylabel("Number of ${O_2}$")
#ax2.set_ylim([100,310])
#ax2.yaxis.set_ticks(np.arange(50, 210, step=50))

ax.set_title(directory[directory.find('LAMMPS'):]+'\n\n'+'Molecule: D')
savedir = 'python_outputs\\figures\\species_plot'+str(random.randint(10000, 99999999))
fig.savefig(savedir, dpi=300, bbox_inches='tight')
plt.show()

print('--------------------------')
ne.print_runtime(time.time()-start_time)
print('--------------------------')