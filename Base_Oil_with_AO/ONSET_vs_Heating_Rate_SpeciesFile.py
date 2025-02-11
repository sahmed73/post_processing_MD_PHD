# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Mar 25 01:06:33 2024
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

def function(x, A, B, C, D):   
    y = C+((D-C)/(1+np.exp((x-A)/B)))
    return y

def inv_function(y,A,B,C,D):
    if D-y<0:
        raise ValueError('Warning from get_onset: Negative D value adjusted!!')
    x = A + B*np.log((D-y)/(y-C))
    return x


baseoil = 'PAO4'
location = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Onset_Ramping_Rate'
# ramps   = ['500Kpps','200Kpps','100Kpps','50Kpps']
ramps   = [200,100,50]

data = []
for ramp in ramps:
    speciesfilepath = location+'\\'+f'{ramp}Kpps'+'\\species.out'
    species = sfp.get_species_count(speciesfilepath)
    data.append(species)
#%%
## Global Variables

formula = {'PAO': 'H62C30', 'Squalane': 'H62C30',
                   'Molecule A':'H30C19O3', 'PAO_Radical': 'H61C30O'}

# Set default font sizes
plt.rcParams['figure.figsize'] = [6, 4]
plt.rc('font', size=12) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=15) # Axes label size
plt.rc('xtick', labelsize=15) # X-axis tick label size
plt.rc('ytick', labelsize=15) # Y-axis tick label size
plt.rc('legend', fontsize=15) # Legend fontsize

pltsave = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\onset'+f'\\onset-{random.randint(0,99999999)}'

timestep     = 0.25
ramp_rate    = 4
initial_temp = 300
skipts       = 0

fig, ax = plt.subplots()
colors = ['tab:red','tab:blue','tab:orange','k','grey','purple','green']
item   = formula['Molecule A']

for i, df in enumerate(data):
    color, ramp = colors[i], ramps[i]
    df = df.T[item]
    df.index = initial_temp + (df.index*timestep*ramp)/(1000)
    print(ramp)
    print(df)
    
    ## creating average sim data out of df
    x = df.index
    y = df.values
    ## curve fit
    popt, cov = curve_fit(function, x, y,p0=[1200,1000,0,25], maxfev = 10000)
    y_fit     = function(x, *popt)
    iimc      = df.iloc[0:5].mean()
    print(f'iimc:{iimc}')
    onset     = np.round(inv_function(iimc-0.5, *popt))    
    ax.plot(df, c=color, label=f'{ramp} K/ps')
    R_squared = mplt.R_square(df, y_fit)
    print(R_squared)
    ax.plot(x,y_fit,linewidth=2,c=color)
    ax.axvline(onset, ymax=0.94, linestyle='--', color=color)
    ax.text(2800,37+3.5*i,f'{np.round(onset)} K', color=color)
    
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('# of intact molecule/radical')
ax.set_ylim(-1,51)
# ax.set_xlim(400,2000)
ax.legend(loc='center left')
ax.grid()
fig.savefig(pltsave,dpi=300,bbox_inches='tight')
