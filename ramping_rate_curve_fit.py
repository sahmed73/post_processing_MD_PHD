# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 06:16:59 2023

@author: arup2
"""

import magnolia.log_parser_FHB as lfp
import magnolia.bondfile_parser as bfp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import random
import sys

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation\AO-1\50_AO1_100_O2\CHO_ReaxFF-2016\Onset\0.5K-per-ps\Sim-1'
filename = '\\log.lammps'
logfile = directory+filename

thermo = lfp.thermo_dict(logfile, 1)
timestep = 0.25
print('Timestep:',timestep)
steps = thermo['Step']
ps = np.array(bfp.step2picosecond(steps, timestep))
temp = thermo['Temp']
energy = np.array(thermo['TotEng'])-min(thermo['TotEng'])
density = thermo['Density']
 
start = bfp.get_timeindex(ps, 0)
end = -1#bfp.get_timeindex(ps, 200)


x = np.array(ps)
y = np.array(temp)
plt.plot(x[start:end],y[start:end])

parameters = lfp.tempramp(thermo, timestep)

m,c = parameters
y_fit = m*x+c

print('Ramping Rate (m):',round(parameters[0],2))
print('Initial Temperature (c):',round(parameters[1],2))
plt.plot(x[start:end],y_fit[start:end],color='r')