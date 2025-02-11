# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Aug 30 09:50:08 2024
"""

############################################
##
#   Heatmap of most abundance species
#   Abundance based on mean number of species
#   at the last 100 ps.
#   Basically, it show only the stable products not intermediates
#   Using species file
##
############################################



import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import random
import pandas as pd
import magnolia.speciesfile_parser as sfp
import seaborn as sns
import numpy as np
import magnolia.plot_template as mplt
import magnolia.access_ucm_cluster as ucm
# mplt.custom_plot_features()


dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\7_A0001-5_eq_NPT-NPT-NVT_Pro_2ramp\A0003\Production\4KPS_Ramp_Rate'

sim_dirr = [''] #,'\\Sim-2','\\Sim-3']
flag = True
cutoff = 0.30
for sim in sim_dirr:
    directory = dirr+sim
    speciesfile_path = directory+'/species.out'
    
    current = sfp.get_species_count(speciesfile_path).T
    if flag:
        output_species = current
        flag = False
    else:
        output_species = output_species.add(current,fill_value=0)

#%%
# User Input
timestep = 0.25
skipts   = 0#(1000-300)/4
nspecies = 10
tailps   = 100  #ps
exclude  = [] #['O2','H30C19O3']
n_sim    = len(sim_dirr)
n_xticks = 6  # how many xtick division in x axis i.e. 1+n_xticks ticklabels 
norm = None#LogNorm(vmin=1, vmax=60)

details  = ('\n\nThis heatmap is generated using species.out file.\nNumber of '
            'species are average over {} simulations.\n'
            'Species are sorted by mean species count of last {} ps.\n'
            'First {} species are showing here.\n'.format(n_sim,tailps, nspecies))

title    = dirr[dirr.find('LAMMPS')+7:]+ details+f"(cutoff: {cutoff})\n\n"



species = output_species.copy()/n_sim # copy is important
# species.replace(0,1e-20,inplace=True)
species.columns = species.columns*timestep/1000-skipts
species    = species.loc[:,0:] # skip the ramping steps
species    = species.drop(exclude)

## mean of last {tailps} ps
# means          = species.loc[:,species.columns.max()-tailps:].mean(axis=1)
# sorted_index   = means.sort_values(ascending=False).index.tolist()
# species        = species.loc[sorted_index]

species       = species.head(nspecies)

# index naming for y axis
species.index = bfp.make_molecular_formula_latex(species.index,sort=True)

new_index = list(species.index)
for i in range(len(new_index)):
    index=new_index[i]
    if index=='$C_{19}H_{30}O_{3}$':
        new_index[i]='(Antioxidant) '+index
    elif index=='$C_{30}H_{61}O$':
        new_index[i]='(PAO Radical) '+index
        
    print(index)
species.index=new_index

xticks_freq = int(species.columns.size/n_xticks)
plt.rcParams['font.size']=14
fig, ax = plt.subplots(figsize=(6,0.38*nspecies),dpi=350)
# species.columns=species.columns*4+600
sns.heatmap(species,cmap='jet',xticklabels=xticks_freq,
                 ax=ax,norm=norm, vmax=20)

# make the xticks integers
xticklabels = ax.get_xticklabels()
xticks = np.array([float(label.get_text()) for label in xticklabels])
ax.set_xticklabels(xticks.astype(int))

ax.set_xlabel('Time (ps)')
plt.tick_params(left=False)
plt.title(title, fontsize=6)
plt.yticks(rotation=0)
plt.grid(False)
# plt.tight_layout()
mplt.saveplot(folder='heatmap_of_species',name='hm_species')
