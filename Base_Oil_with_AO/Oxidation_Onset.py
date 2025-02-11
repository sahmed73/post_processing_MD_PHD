# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Feb 25 15:23:13 2024
"""
import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import os
import random
import sys
import pandas as pd
import magnolia.plot_template as mplt
from scipy.optimize import curve_fit

baseoil = 'PAO4'
location = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Onset'
sim_dir = ['Sim-1','Sim-2','Sim-3']    

atomsymbols= ['H','C','O']
bonddata = {}
for sim in sim_dir:
    bondfilepath = location+'\\'+sim+'\\bonds.reaxc'
    bonddata[sim]= bfp.parsebondfile(bondfilepath,cutoff=0.3)
#%%
## Global Variables
timestep        = 0.25
ramp_rate       = 4
initial_temp    = 300
baseoil_formula = {'PAO4': 'H62C30', 'Squalane': 'H62C30'}

df = pd.DataFrame()
for sim in sim_dir:
    neighbours = bonddata[sim]['neighbours']
    atypes    = bonddata[sim]['atypes']
    
    CO_bond_count = {}
    for step, neigh in neighbours.items():
        count = 0
        ps = initial_temp + step/(timestep*ramp_rate*1000)
        for parent, children in neigh.items():
            # if parent is oxygen and childrend have carbon
            if atypes[parent]==3:
                childrend_type = [atypes[u] for u in children]
                count+=childrend_type.count(2)
        CO_bond_count[ps]=count
    df[sim]=pd.Series(CO_bond_count)
#%%
idx = df.mean(axis=1).idxmax()
bdf = df
# Set default font sizes
plt.rc('font', size=15) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=16) # Axes label size
plt.rc('xtick', labelsize=15) # X-axis tick label size
plt.rc('ytick', labelsize=15) # Y-axis tick label size
plt.rc('legend', fontsize=12) # Legend fontsize

pltsave = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\oxidation'+f'\\oxidation-{random.randint(0,99999999)}'
title = location[location.find('Base'):]

fig, ax = plt.subplots()
mplt.mean_range_plot(bdf,color='blue',ax=ax)
ax.set_xlabel('Temperature')
ax.set_ylabel('Number of C-O bonds')
ax.set_title(title)

def gaussian(x, a, b, c, d):
    return a * np.exp(-0.5 * ((x - b) / c) ** 2) + d

temp = df.index
guess = (df.mean(axis=1).max(), df.mean(axis=1).idxmax(), 100, 25)
popt, pcov = curve_fit(gaussian, temp, df.mean(axis=1), p0=guess)
y_fit = gaussian(temp,*popt)

R_squared = mplt.R_square(df.mean(axis=1), y_fit)
print(R_squared)

y_initial = y_fit[0]
for y,t in zip(y_fit,temp):
    if abs(y_initial-y)>=2:
        print(y,t)
        oxidation_onset = t
        break

# ax.plot(temp,gaussian(temp,*popt),color='r')
ax.set_xlim(1250,idx)
ax.axvline(oxidation_onset, ymax=0.5, linestyle='--',color='b')
ax.text(oxidation_onset-100, 135, f'{oxidation_onset} K', color='b')
fig.savefig(pltsave,bbox_inches='tight',dpi=300)