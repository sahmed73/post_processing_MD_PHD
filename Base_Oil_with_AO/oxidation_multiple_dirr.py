# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Feb 26 09:06:13 2024
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
locations = [r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\0.1fs_25_PAO4_80_O2_Kowalik_2019\Onset',
             r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\0.1fs_25_PAO4_80_O2_6A_Kowalik_2019\Onset']
##################-important-###############
N = [25,25]
data = []
for i, location in enumerate(locations):
    sim_dir = ['Sim-1']#,'Sim-2','Sim-3']    
    atomsymbols= ['H','C','O']
    bonddata = {}
    for sim in sim_dir:
        bondfilepath = location+'\\'+sim+'\\bonds.reaxc'
        bonddata[sim]= bfp.parsebondfile(bondfilepath,cutoff=0.3,mtypes=True)
    
    df = pd.DataFrame()
    for sim in sim_dir:
        neighbours = bonddata[sim]['neighbours']
        atypes     = bonddata[sim]['atypes']
        mtypes     = bonddata[sim]['mtypes']
        
        CO_bond_count = {}
        for step, neigh in neighbours.items():
            count = 0
            # ps   = (step*timestep)/1000
            # temp = initial_temp + ps*ramp_rate
            for parent, children in neigh.items():
                # if parent is oxygen and childrend have carbon
                if atypes[parent]==3:
                    for child in children:
                        if 1<=mtypes[child]<=N[i] and atypes[child]==2:
                            count+=1
            CO_bond_count[step]=count
        df[sim]=pd.Series(CO_bond_count)
    data.append(df)
#%%
timestep        = 0.1
ramp_rate       = 4
initial_temp    = 300
baseoil_formula = {'PAO4': 'H62C30', 'Squalane': 'H62C30'}

# Set default font sizes
plt.rcParams['figure.figsize'] = [6, 4]
plt.rc('font', size=14) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=17) # Axes label size
plt.rc('xtick', labelsize=16) # X-axis tick label size
plt.rc('ytick', labelsize=16) # Y-axis tick label size
plt.rc('legend', fontsize=12) # Legend fontsize
pltsave = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\oxidation'+f'\\oxidation-{random.randint(0,99999999)}'

def gaussian(x, a, b, c, d):
    return a * np.exp(-0.5 * ((x - b) / c) ** 2) + d

skipts = (1600-300)/ramp_rate
colors = ['tab:blue','tab:red','tab:orange']
labels = ['25 PAO+80 O$_2$', '25 PAO+6 A+80 O$_2$']

fig, ax = plt.subplots()
for i, df in enumerate(data):
    idx = df.mean(axis=1).idxmax() 
    adf = df.copy()
    adf.index = initial_temp+df.index*timestep*ramp_rate/1000
    bdf = adf.loc[0:,:].copy()
    temp= bdf.index
    
    mplt.mean_range_plot(bdf,ax=ax, label=labels[i], color=colors[i])
    # guess = (bdf.mean(axis=1).max(), bdf.mean(axis=1).idxmax(), 100, 25)
    # # popt, pcov = curve_fit(gaussian, temp, bdf.mean(axis=1), p0=guess)
    # # y_fit = gaussian(temp,*popt)
    
    # R_squared = mplt.R_square(bdf.mean(axis=1), y_fit)
    # print(R_squared)
    
    # y_initial = y_fit[0]
    # for y,t in zip(y_fit,temp):
    #     if abs(y_initial-y)>=2:
    #         print(y,t)
    #         oxidation_onset = t
    #         break
    
    # ax.plot(temp,gaussian(temp,*popt))
    # # ax.set_xlim(0,400)
    # ax.axvline(oxidation_onset, ymax=0.5, linestyle='--', color=colors[i])
    # # ax.text(oxidation_onset-100, 135, f'{oxidation_onset} K', color=colors[i])

ax.set_xlim(1500,2500)
ax.set_ylim(-1,50)
ax.set_xticks(np.arange(1500,2501,200))
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Number of C-O bonds')
ax.legend(loc='upper left')
fig.savefig(pltsave,bbox_inches='tight',dpi=300)