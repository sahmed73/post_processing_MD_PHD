# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 13 23:14:09 2023
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import magnolia.plot_template as mplt
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import networkx as nx

fig, ax = plt.subplots()
double = 1.2 # double bond
unsaturate = pd.DataFrame()
timestep = 0.25
initial_temp  = 300
temp_ramprate = 4
colors = ['gold', 'green']
oils = ['PAO4','Squalane']
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
        atomConnectivity = bfp.parsebondfile_asGraph((bondfilepath))
        ring_count = []
        for step, atomConnectGraphs in atomConnectivity.items():
            cycles = nx.cycle_basis(atomConnectGraphs)
            benzene = [x for x in cycles if len(x)==6]
            ring_count.append(len(benzene))
        
        df[sim]=ring_count
        
    # df.index=(steps*timestep/1000)*temp_ramprate+initial_temp
    data.append(df)

#%% ploting using own function
for i, df in enumerate(data):
    df.index=range(301,3501)
    dff = df.loc[:3300,:]
    mplt.mean_range_plot(dff, color=colors[i],label=oils[i],
                         alpha=0.0,linewidth=0.5)
    plt.legend(loc='upper left')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Number of benzene cycle')
  
title = ",".join(oils+simulations)+'\n'
plt.title(title,fontsize=10)
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures'
plt.savefig(savedir+'\\benzene_count{}'.format(random.randint(0, 10000000)), dpi=500)
plt.show()