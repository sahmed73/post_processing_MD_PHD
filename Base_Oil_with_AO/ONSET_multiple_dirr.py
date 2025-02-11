# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Feb 26 09:16:29 2024
"""

import magnolia.bondfile_parser as bfp
import magnolia.speciesfile_parser as sfp
import numpy as np
import matplotlib.pyplot as plt
import os
import random
import sys
import pandas as pd
import magnolia.plot_template as mplt
from scipy.optimize import curve_fit

baseoil = 'PAO4'
locations = [r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Onset_Ramping_Rate']
temps   = ['500Kpps','200Kpps','100Kpps','50Kpps']

data = []
for location in locations:
    sim_dir = ['']#['Sim-1']#,'Sim-2','Sim-3']
    atomsymbols= ['H','C','O']
    species = {}
    for sim in sim_dir:
        bondfilepath    = location+'\\'+sim+'\\bonds.reaxc'
        speciesfilepath = location+'\\'+sim+'\\species.out'
        species[sim] = sfp.get_species_count(speciesfilepath)
        # species[sim] = bfp.get_species_count(bondfilepath, atomsymbols)
    data.append(species)
#%%
## Global Variables

baseoil_formula = {'PAO': 'H62C30', 'Squalane': 'H62C30',
                   'Molecule A':'H30C19O3', 'PAO_Radical': 'H61C30O'}

# Set default font sizes
plt.rcParams['figure.figsize'] = [6, 4]
plt.rc('font', size=18) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=15) # Axes label size
plt.rc('xtick', labelsize=15) # X-axis tick label size
plt.rc('ytick', labelsize=15) # Y-axis tick label size
plt.rc('legend', fontsize=15) # Legend fontsize

pltsave = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\onset'+f'\\onset-{random.randint(0,99999999)}'

def function(x, A, B, C, D):   
    y = C+((D-C)/(1+np.exp((x-A)/B)))
    return y

def inv_function(y,A,B,C,D):
    if D-y<0:
        raise ValueError('Warning from get_onset: Negative D value adjusted!!')
    x = A + B*np.log((D-y)/(y-C))
    return x

fig, ax = plt.subplots()
colors = ['tab:red','tab:blue','tab:orange']
labels = ['PAO Radical','PAO']

oils = ['PAO_Radical', 'PAO', 'PAO_Radical']
for i, species in enumerate(data):
    timestep = 0.25
    skipts   = 0
    
    if i == 1: 
        timestep=0.10
    ramp_rate       = 4
    initial_temp    = 300
    
    df = pd.DataFrame()
    for sim in sim_dir:
        spec       = species[sim].loc[baseoil_formula[oils[i]]]
        spec.index = initial_temp + (spec.index*timestep*ramp_rate)/(1000)
        df[sim]    = spec
    
    ## creating average sim data out of df
    x = df.index
    y = df.mean(axis=1)
    ## curve fit
    popt, cov = curve_fit(function, x, y,p0=[2200,1000,0,25], maxfev = 10000)
    y_fit     = function(x, *popt)
    iimc      = df.iloc[0,:].mean()
    onset     = np.round(inv_function(iimc-1, *popt))
    error = np.round(df.mean(axis=1).std())
    print(f'error {error}%')   
    
    color, label = colors[i], labels[i]
    mplt.mean_range_plot(df, ax=ax, meanplot=True, shade=False, c=color,
                         label=label)
    R_squared = mplt.R_square(df.mean(axis=1), y_fit)
    print(R_squared)
    ax.plot(x,y_fit,linewidth=2,c=color)
    ax.axvline(onset, ymax=0.94, linestyle='--', color=color)
    ax.text(onset-150,5,f'{np.round(onset)} K', rotation=90, color=color)
    
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('# of intact molecule/radical')
ax.set_ylim(-1,26)
# ax.set_xlim(400,2000)
ax.legend()
ax.grid()
fig.savefig(pltsave,dpi=300,bbox_inches='tight')
