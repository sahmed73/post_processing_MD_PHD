# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 16 11:15:52 2023
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd
import sys
from scipy.optimize import curve_fit, newton

oil = 'Squalane'

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


path = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Onset'.format(oil,oil)
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
df2 = pd.DataFrame()
for sim in ["Sim-1","Sim-2", "Sim-3"]:
### Directory ###
    directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Onset\{}".format(oil,oil,sim)
    filename  = "\\bonds.reaxc"
    bondfilepath = directory+filename
    
    ### Parsing Bondfile ###
    bonddata = bfp.parsebondfile(bondfilepath,mtypes=True)
    
    ## Geting neighbour lists
    neighbours = bonddata['neighbours']
    atypes     = bonddata['atypes']
    mtypes     = bonddata['mtypes']
    asyms      = ['H','C','O']
    steps      = np.array(list(neighbours.keys()))
    ### Loop through the neighbours ###
    COBond_count = []
    for step, neigh in neighbours.items():
        count = 0
        for parent, children in neigh.items():
            if atypes[parent]==2:
                count += [atypes[child] for child in children].count(3)
        COBond_count.append(count)
    
    df2[sim]=COBond_count
#%%
color= 'green'
## Getting the Upper and Lower bound
upper_bound = df.iloc[:, -3:].max(axis=1)
lower_bound = df.iloc[:, -3:].min(axis=1)

## plotting the fill-between plot
x       = df['temp']
fig, ax = plt.subplots()
ax.fill_between(x, lower_bound, upper_bound, color=color,alpha=0.2)

## creating average sim data out of df
y         = df.iloc[:,-3:].mean(axis=1)

## curve fit
popt, cov = curve_fit(function, x, y,p0=[2200,107,0,25])
y_fit     = function(x, *popt)

## plot fit values
ax.plot(x,y_fit,color=color,label='Degradation')

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
# ax.plot([onset]*2,[-0.5,y_target],'--',color=color)
                
df2['mean'] = df2.iloc[:,-3:].mean(axis=1)
print(df2)


ax.scatter(x,df2['mean'],label='Oxidation',color=color,s=1)

## plot
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Number of moleculs')
ax.set_ylim(-5)
ax.set_xlim(300,3500+100)
plt.legend(loc='upper left')
fig.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\onset\onset_{}'.format(random.randint(0,10000000)),dpi=400)