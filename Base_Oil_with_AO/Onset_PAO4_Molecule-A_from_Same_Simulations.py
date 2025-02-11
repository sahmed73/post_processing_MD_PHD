# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Feb 26 09:16:29 2024
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
locations = [r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_6A_Soria\Production\1600']

data = []
for location in locations:
    sim_dir = ['Sim-1','Sim-2','Sim-3']
    atomsymbols= ['H','C','O']
    species = {}
    for sim in sim_dir:
        bondfilepath = location+'\\'+sim+'\\bonds.reaxc'
        species[sim] = bfp.get_species_count(bondfilepath, atomsymbols)
    data.append(species)
#%%
## Global Variables
timestep        = 0.25
ramp_rate       = 4
initial_temp    = 300
baseoil_formula = {'PAO4': 'H62C30', 'Squalane': 'H62C30',
                   'Molecule A':'H30C19O3'}

# Set default font sizes
plt.rc('font', size=12) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=15) # Axes label size
plt.rc('xtick', labelsize=15) # X-axis tick label size
plt.rc('ytick', labelsize=15) # Y-axis tick label size
plt.rc('legend', fontsize=12) # Legend fontsize
pltsave = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\onset'+f'\\onset-{random.randint(0,99999999)}'

def function(x, A, B, C, D):   
    y = C+((D-C)/(1+np.exp((x-A)/B)))
    return y

def inv_function(y,A,B,C,D):
    if D-y<0:
        raise ValueError('Warning: Negative D value encounterd!')
    x = A + B*np.log((D-y)/(y-C))
    return x

fig, ax = plt.subplots()
colors = ['tab:blue','tab:orange','tab:red']
labels = ['PAO 4', 'Molecule A']

oils = ['PAO4','Molecule A']

species = data[0]
vlines  = [305,250]
for i, oil in enumerate(oils):
    df = pd.DataFrame()
    for sim in sim_dir:
        spec       = species[sim].loc[baseoil_formula[oil]].copy()
        aspec      = spec.copy()
        aspec.index= spec.index*timestep/1000
        bspec      = aspec.loc[0:].copy()
        df[sim]    = bspec.copy()
    
    ## creating average sim data out of df
    x = df.index
    y = df.mean(axis=1)
    ## curve fit
    popt, cov = curve_fit(function, x, y,p0=[2200,107,0,25])
    y_fit     = function(x, *popt)
    iimc      = df.iloc[0,:].mean()
    onset     = np.round(inv_function(iimc-1, *popt))
    error = np.round(df.mean(axis=1).std())
    print(f'error {error}%')   
    
    mplt.mean_range_plot(df, ax=ax, color=colors[i], label=labels[i])
    R_squared = mplt.R_square(df.mean(axis=1), y_fit)
    print(R_squared)
    # ax.set_ylim(-1,26)
    # ax.plot(x,y_fit,linewidth=2,color=colors[i],label=labels[i])
    ax.axvline(vlines[i], ymax=0.94, linestyle='--',color=colors[i])
    # ax.text(onset-90,10+6*i,f'{np.round(onset)} K', rotation=90, color=colors[i])
    
ax.set_xlabel('Temperature')
ax.set_ylabel('Number of intact molecule')
ax.set_ylim(-1,26)
ax.legend()
fig.savefig(pltsave,dpi=300,bbox_inches='tight')
