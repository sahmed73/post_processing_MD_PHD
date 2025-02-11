# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Feb 25 21:56:25 2024
"""

import magnolia.speciesfile_parser as sfp
import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
import random
import pandas as pd
import numpy as np
import magnolia.plot_template as mplt
import matplotlib as mpl
plt.style.use('classic')
plt_width = 7
aspect_ratio = 1.333333
plt.figure(figsize=[plt_width,plt_width/aspect_ratio])
plt.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.labelpad'] = 15
plt.rc('axes', grid=True)
plt.rc('font', size=10) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=20) # Axes label size
plt.rc('xtick', labelsize=20) # X-axis tick label size
plt.rc('ytick', labelsize=20) # Y-axis tick label size
plt.rc('legend', fontsize=20) # Legend fontsize


base_oil = 'PAO4'
    
location = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO_Radical\20_PAO_Radical_20_A_Soria\Production\800K"

sim_dir = ['Sim-1'] # ['Sim-1', 'Sim-2','Sim-3']
species = {}
for sim in sim_dir:
    speciesfilepath = location+"\\"+sim+'\\species.out'
    species[sim] = sfp.get_species_count(speciesfilepath)


saveplt = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\species'+f'\\timeseries_species_{random.randint(0,100000000)}'

timestep        = 0.25
ramp_rate       = 4
initial_temp    = 300
skipts          = 0#(900-initial_temp)/ramp_rate

order   =  ['H30C19O3']#['H61C30O', 'H24C15O',]# ['H24C15O', 'H61C30O']#
order_latex = bfp.make_molecular_formula_latex(order,sort=True)
mapping =  dict(zip(order,order_latex))


for key in ['H30C19O3', 'H61C30O', 'H24C15O', 'H59C29']:
    if key in mapping:
        if key=='H30C19O3': mapping[key] +='(A)'
        if key=='H61C30O': mapping[key]  += '(R-CH2O)'#'(PAO radical)'
        if key=='H24C15O': mapping[key]  += '(BHT)'
        if key=='H59C29': mapping[key]  += '(R-)'

species = species[sim_dir[0]]
species = species.sort_values(by=species.columns[-1], ascending=False)
df = species.T
df = df[order]
print(df.columns)
df.index = df.index*timestep/1000-skipts
df = df.iloc[:,:10]

plt.figure(figsize=[6,4])
# labels = bfp.make_molecular_formula_latex(df.columns,sort=True)
# labels = [mapping[x] for x in df.columns]
plt.plot(df,label=df.columns[0])
plt.legend(loc=[1,0])
plt.xlabel('Time (ps)')
plt.ylabel('Number of molecule')
plt.grid()
# plt.xlim(0,18)
plt.ylim(-1,21)
plt.savefig(saveplt,dpi=300,bbox_inches='tight')