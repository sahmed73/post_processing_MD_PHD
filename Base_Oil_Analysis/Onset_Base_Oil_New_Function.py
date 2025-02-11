# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 16 01:03:47 2023
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd
import sys
from scipy.optimize import curve_fit, newton


## plotting fit curve
def function(x, A, B, C, D):   
    y = C+((D-C)/(1+np.exp((x-A)/B)))
    return y

def inv_function(y,A,B,C,D):
    if D-y<0:
        print('Warning from get_onset: Negative D value adjusted!!')
        D = df.iloc[:100,-1].mean()
    x = A + B*np.log((D-y)/(y-C))
    return x


path = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Squalane\25_Squalane_200_O2_Soria\Onset'
whole = 'H62C30'
atomsymbols = ['H','C','O']
timestep = 0.25 # fs
temp_ramp = 4 # K/ps
initial_temp = 300

# fig, ax = bfp.onset_plot(path, whole, atomsymbols,
#                          timestep, temp_ramp, initial_temp, color='r')


df         = pd.DataFrame()
onset_list = []
## Iterate over number of simulations
for sim in ['Sim-1','Sim-2','Sim-3']:
    bondfilepath = path+'\\'+sim+'\\bonds.reaxc'

    ## Geting neighbours for each simulations
    bonddata   = bfp.parsebondfile(bondfilepath)
    neighbours = bonddata['neighbours']
    atypes     = bonddata['atypes']
    
    ## Getting temperatures using steps, timesteps, temp_ramp, initial_temp
    steps      = np.array(list(neighbours.keys()))
    time       = steps*timestep/1000  # in piccosecond
    temp       = initial_temp + time*temp_ramp
    df['temp'] = temp
    
    ##  Getting number of whole
    nwhole     = []
    for step, neigh in neighbours.items():
        count = 0
        molecules = bfp.get_molecules(neigh)
        for molecule in molecules:
            species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
            if species == whole:
                count+=1
        nwhole.append(count)
        
    nwhole = np.array(nwhole)
    df[sim]=nwhole
    
    ## getting individual onset
    popt, _    = curve_fit(function,temp,nwhole, p0=[2200,107,0,25])
    nwhole_fit = function(temp, *popt)
    onset      = inv_function(24, *popt)
    onset_list.append(onset)
    
    ## finding r square value
    nwhole_mean = np.mean(nwhole)
    tss = np.sum((nwhole - nwhole_mean)**2)
    rss = np.sum((nwhole - nwhole_fit)**2)
    r2 = 1 - (rss / tss)
    
    ## printing sim data
    print(sim)
    print('Indivisual onset:',onset)
    print("R squared value:",r2)
    print()

#%%
color= 'green'
## Getting the Upper and Lower bound
upper_bound = df.iloc[:, -3:].max(axis=1)
lower_bound = df.iloc[:, -3:].min(axis=1)

## plotting the fill-between plot
x       = df['temp']
fig, ax = plt.subplots()
ax.fill_between(x, lower_bound, upper_bound, color=color,alpha=0.2)
ax.set_ylim(-0.5,27)



## creating average sim data out of df
y         = df.iloc[:,-3:].mean(axis=1)

## curve fit
popt, cov = curve_fit(function, x, y,p0=[2200,107,0,25])
y_fit     = function(x, *popt)

## plot fit values
ax.plot(x,y_fit,color=color,label='Squalane')

## finding r square value
y_mean = np.mean(y)
tss = np.sum((y - y_mean)**2)
rss = np.sum((y - y_fit)**2)
r2 = 1 - (rss / tss)

## getting onset 
y_target  = df.iloc[:100,-1].mean()-1
onset     = inv_function(y_target, *popt)

## geting std/error
std = np.array(nwhole).std()

print('Overall Onset: {} ± {} K'.format(onset,std))
print('R squared value:',r2)

## plotting a onset indicating verticle dashed line
ax.plot([onset]*2,[-0.5,y_target],'--',color=color)

## plot
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Number of moleculs')
plt.legend()
fig.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\onset\onset_{}'.format(random.randint(0,10000000)),dpi=400)