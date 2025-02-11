# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Feb 23 10:13:35 2024
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

baseoil = 'PAO_Radical'
location = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\Pure_50B\Onset'
sim_dir = ['Sim-1']#,'Sim-2','Sim-3']    

atomsymbols= ['H','C','O']
species = {}
for sim in sim_dir:
    bondfilepath    = location+'\\'+sim+'\\bonds.reaxc'
    speciesfilepath = location+'\\'+sim+'\\species.out'
    species[sim] = sfp.get_species_count(speciesfilepath).T
    # species[sim] = bfp.get_species_count(bondfilepath, atomsymbols)
#%%
## Global Variables
timestep        = 0.25
ramp_rate       = 4
initial_temp    = 300


df = pd.DataFrame()
for sim in sim_dir:
    spec       = species[sim].loc['H34C26O4']
    spec.index = initial_temp + (spec.index)/(timestep*ramp_rate*1000)
    df[sim]    = spec
    
def function(x, A, B, C, D):   
    y = C+((D-C)/(1+np.exp((x-A)/B)))
    return y

def inv_function(y,A,B,C,D):
    if D-y<0:
        raise ValueError('Warning from get_onset: Negative D value adjusted!!')
    x = A + B*np.log((D-y)/(y-C))
    return x

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
title = location[location.find('Base'):]+'\n'+f'onset {onset}$\pm${error}\n'

# Set default font sizes
plt.rc('font', size=22) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=22) # Axes label size
plt.rc('xtick', labelsize=20) # X-axis tick label size
plt.rc('ytick', labelsize=20) # Y-axis tick label size
plt.rc('legend', fontsize=12) # Legend fontsize

color = 'tab:red'
fig, ax = plt.subplots()
mplt.mean_range_plot(df, ax=ax, meanplot=True, shade=False, c=color)
R_squared = mplt.R_square(df.mean(axis=1), y_fit)
print(R_squared)
# ax.set_ylim(-1,21)
# ax.set_xlim(400,2000)
ax.set_title(title)
ax.plot(x,y_fit,linewidth=2,c=color)
ax.axvline(onset, ymax=0.94, linestyle='--', color=color)
# ax.text(onset-150,5,f'{np.round(onset)} K', rotation=90, color=color)
ax.set_xlabel('Temperature')
ax.set_ylabel('Number of intact molecule')
xticks = range(500,2001,250)
# ax.set_xticks(xticks)
ax.grid()
mplt.saveplot('onset','Pure_B_Onset')
