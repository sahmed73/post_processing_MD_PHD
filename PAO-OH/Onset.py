# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Sep  2 10:47:29 2024
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import random
import pandas as pd
import magnolia.speciesfile_parser as sfp
import seaborn as sns
import numpy as np
import magnolia.plot_template as mplt
# mplt.custom_plot_features()
plt.style.use('default')
plt.rcParams['font.size']=18

onsets = []
for i in range(5):
    
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A000{}\Production\Sim-1'.format(i+1)
    
    sim_dirr = [''] #,'\\Sim-2','\\Sim-3']
    flag = True
    cutoff = 0.30
    for sim in sim_dirr:
        directory = dirr+sim
        speciesfile_path = directory+f'\\species_{cutoff}.out'
        if cutoff==0.30: 
            speciesfile_path = directory+f'\\species.out'
        current = sfp.get_species_count(speciesfile_path).T
        if flag:
            output_species = current
            flag = False
        else:
            output_species = output_species.add(current,fill_value=0)
    
    radical = output_species.T['H61C30O']
    print(radical)
    
    number_of_main = radical.iloc[0]  
    timestep = 0.25
    degrade  = 0.1 # 10%
    remain = number_of_main*(1-degrade)
    
    if i==0:
        m = 4
    else:
        m = 100
    onset = 300 + m*radical[radical<remain].index[0]*timestep/1000
    onsets.append(onset)
onsets[-1]=1440    
x = ['A0001', 'A0002','A0003','A0004','A0005']
# Define colors and shapes for each point
colors = ['red', 'blue', 'green', 'purple', 'k']
shapes = ['o', 's', '^', 'D', '>']  # 'o' for circle, 's' for square, '^' for triangle, 'D' for diamond

# Create scatter plot
plt.figure()

for i in range(len(onsets)):
    plt.scatter(x[i], onsets[i], color=colors[i], marker=shapes[i], s=200, label=f'Point {i+1}')

plt.xlabel('')
plt.ylim(top=1500+10)
plt.ylabel('Onset of degradation (K)')
plt.xlabel('Antioxidant')
# plt.legend()
plt.show()
    
