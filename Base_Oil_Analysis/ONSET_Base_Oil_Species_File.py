# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 13 13:14:24 2023
"""

import magnolia.bondfile_parser as bfp
import magnolia.speciesfile_parser as sfp
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from scipy.optimize import curve_fit


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

fig, ax = plt.subplots()
colors= ['b','g']

for i, baseoil in enumerate(['PAO4','Squalane']):
    path = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Onset'.format(baseoil,baseoil)
    whole = 'H62C30'
    atomsymbols = ['H','C','O']
    timestep = 0.25 # fs
    temp_ramp = 4 # K/ps
    initial_nmolecule = 25
    initial_temp = 300
    
    onset_list = []
    df = pd.DataFrame()
    ## Iterate over number of simulations
    for sim in ['Sim-1','Sim-2','Sim-3']:
        speciesfilepath = path+'\\'+sim+'\\species.out'
    
        ## Geting number of species dataframe for each simulations
        nspecies = sfp.get_species(speciesfilepath).T
        nspecies.index = nspecies.index * 0.25*4/1000 + 300
        df[sim]=nspecies[whole]
        
        ## getting individual onset
        popt, _    = curve_fit(function,df.index,nspecies[whole], p0=[2200,107,0,25])
        nwhole_fit = function(df.index, *popt)
        onset      = round(inv_function(24, *popt),0)
        onset_list.append(onset)
        
        ## finding r square value
        tss = np.sum((nspecies[whole] - nspecies[whole].mean())**2)
        rss = np.sum((nspecies[whole] - nwhole_fit)**2)
        r2 = round(1 - (rss / tss),4)
        
        ## printing sim data
        print(sim)
        print('Indivisual onset:',onset)
        print("R squared value:",r2)
        print()
    
    ## Getting the Upper and Lower bound
    upper_bound = df.iloc[:, -3:].max(axis=1)
    lower_bound = df.iloc[:, -3:].min(axis=1)
    
    ## plotting the fill-between plot
    x       = df.index
    ax.fill_between(x, lower_bound, upper_bound, color=colors[i],alpha=0.2)
    ax.set_ylim(-0.5,initial_nmolecule+2)
    
    
    
    ## creating average sim data out of df
    y         = df.iloc[:,-3:].mean(axis=1)
    
    ## curve fit
    popt, cov = curve_fit(function, x, y,p0=[2200,107,0,25])
    y_fit     = function(x, *popt)
    
    ## plot fit values
    ax.plot(x,y_fit,color=colors[i],label=baseoil)
    
    ## finding r square value
    y_mean = np.mean(y)
    tss = np.sum((y - y_mean)**2)
    rss = np.sum((y - y_fit)**2)
    r2 = round(1 - (rss / tss),4)
    
    ## getting onset 
    y_target  = df.iloc[:100,-1].mean()-1
    onset     = round(inv_function(y_target, *popt),0)
    
    ## geting std/error
    std = round(np.array(nspecies[whole]).std(),0)
    
    print('Overall Onset: {} ± {} K'.format(onset,std))
    print('R squared value:',r2)
    
    ## plotting a onset indicating verticle dashed line
    ax.plot([onset]*2,[-0.5,y_target],'--',color=colors[i])

## plot
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Number of moleculs')
plt.legend()
fig.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\onset\onset_{}'.format(random.randint(0,10000000)),dpi=400)