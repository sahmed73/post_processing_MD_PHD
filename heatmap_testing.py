# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 13:51:11 2023

@author: Shihab
"""

import magnolia.bondfile_parser as bfp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import random
import sys
#%%-------------------
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO2\50_AO2_300_O2\Production\Sim-1'

filename  = '\\bonds.reaxc'
bondfile  = directory+filename

timestep = 0.25
neighbours,atomtypes = bfp.get_neighbours(bondfile)
steps = list(neighbours.keys())
atomsymbols = ['H','C','O'] # same order as atom type

data = bfp.get_SpeciesCountAtEveryTimestep(neighbours, atomtypes, 'HCO',step2ps=True,step2psargs=[0.25])
#data = bfp.stepwise_species_count(neighbours, atomtypes, atomsymbols,step2ps=timestep)

#%%-------------------------------
number_of_molecules = 25

df = pd.DataFrame(data).fillna(1)
df.index = bfp.make_molecular_formula_latex(df.index,sort=True)

ranked_by_mean_abundance = ''
specie_order = ['$C_{26}H_{34}O_{4}$','$O_{2}$', '$H_{2}O$',  '$CH_{2}O$', '$CO_{2}$', '$CH_{4}O$', '$C_{26}H_{33}O_{4}$', '$C_{26}H_{32}O_{4}$', '$C_{26}H_{33}O_{3}$', '$C_{26}H_{32}O_{3}$', '$H_{2}O_{2}$', '$CH_{2}O_{2}$', '$C_{26}H_{31}O_{3}$', '$CHO_{2}$', '$C_{26}H_{31}O_{4}$', '$CO$', '$C_{26}H_{29}O_{4}$', '$C_{26}H_{30}O_{3}$', '$CH_{4}O_{2}$', '$C_{26}H_{34}O_{5}$', '$C_{26}H_{29}O_{3}$', '$C_{26}H_{30}O_{4}$', '$C_{26}H_{35}O_{5}$', '$C_{26}H_{32}O_{2}$', '$H_{2}$']#['$C_{10}H_{12}O_{3}$','$O_{2}$', '$CO$', '$H_{2}O$', '$CH_{2}O$', '$HO$', '$H_{2}$', '$H$', '$C_{2}H_{2}$', '$C_{10}H_{10}O_{2}$', '$CH_{3}$', '$CO_{2}$', '$HO_{2}$', '$CH_{4}O$', '$C_{2}O$', '$O$', '$C_{10}H_{11}O_{2}$', '$CH_{3}O$', '$C_{10}H_{11}O_{3}$', '$C_{2}H_{2}O$', '$C_{2}$', '$C_{2}HO$', '$C_{9}H_{9}O_{2}$', '$C_{9}H_{8}O_{2}$', '$C_{2}H_{3}$']
if specie_order:
    df = df.loc[specie_order,:]
else: #ranked by mean abundance
    #ranked_by_mean_abundance = 'ranked by mean abundance'
    df.loc[:,'sum'] = df.sum(axis=1)
    df = df.sort_values(by='sum',ascending=False)
    df = df.drop(['sum'],axis=1)

df = df.head(number_of_molecules)
print('Species ORDER:',list(df.index))

_, ax1 = plt.subplots(figsize=(15,8))
ax = sns.heatmap(df,cmap='jet',ax=ax1 ,cbar_kws={'label': 'Number of molecules'},norm=LogNorm(),xticklabels=1000)


ax.set_title(directory[directory.find('AO_Oxidation'):]+'\n\n'+'Top {} molecules '.format(number_of_molecules)+ranked_by_mean_abundance)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Species')
plt.tick_params(left=False,bottom=False)
        

fig = ax.get_figure()
savedir = 'python_outputs\\heatmaps\\heatmap'+str(random.randint(100000, 999999))
fig.savefig(savedir, dpi=300, bbox_inches='tight')
plt.show()

