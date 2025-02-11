# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Oct 12 00:03:48 2023
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd
import sys
from scipy.optimize import curve_fit
import magnolia.plot_template as mplt

df = pd.DataFrame()
double = 1.2 # double bond
triple = 1.9 # triple bond
timestep = 0.25
initial_temp  = 300
temp_ramprate = 4
colors = ['gold', 'green']
oils = ['PAO4']#,'Squalane']
simulations = ["Sim-1","Sim-2", "Sim-3"]

data = []
for oil in oils:
    df = pd.DataFrame()
    for sim in simulations:
        print(oil,sim)
    ### Directory ###
        directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Onset\{}".format(oil,oil,sim)
        filename  = "\\bonds.reaxc"
        bondfilepath = directory+filename
        
        ### Parsing Bondfile ###
        bonddata = bfp.parsebondfile(bondfilepath,bo=True)
        
        ## Geting neighbour lists
        neighbours = bonddata['neighbours']
        atypes     = bonddata['atypes']
        bondorders = bonddata['bondorders']
        asyms      = ['H','C','O']
        steps      = np.array(list(neighbours.keys()))
        ### Loop through the neighbours ###
        COBond_count = []
        for step, neigh in neighbours.items():
            count = 0
            for parent, children in neigh.items():
                if atypes[parent]==2:
                    for child in children:
                        bo = bondorders[step][parent][child]
                        if atypes[child]==3:
                            count+=1
            COBond_count.append(count)
        
        df[sim]=COBond_count
    
    df.index=(steps*timestep/1000)*temp_ramprate+initial_temp
    data.append(df)
    
#%% ploting 
for i, df in enumerate(data):
    mplt.mean_range_plot(df.loc[1500:2250,:], color='b',
                         alpha=0.1,linewidth=0.7)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Number of unsatuareted bonds')

# plt.plot([1623,1623],[0,80],'--')    
title = directory[directory.find('Base_Oil'):-6]+'\n'
plt.xlabel('Temperature (K)')
plt.ylabel('Number of C-O bonds')
plt.title(title)
# plt.xlim(0,3520)
plt.ylim(-1)
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\oxidation'
plt.savefig(savedir+'\\Oxidation_{}'.format(random.randint(0, 10000000)), dpi=400)