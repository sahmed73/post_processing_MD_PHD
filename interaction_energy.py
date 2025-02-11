# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Jan 31 13:39:42 2025
"""

import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import numpy as np

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=18


dirrs = [r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0001\Eq\Sim-3",
         r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0003\Eq\Sim-9_Interaction_Energy"]
filename = r'\log.lammps'
labels = ['A0001', 'A0003']
fig, ax = plt.subplots(dpi=300)
for i, dirr in enumerate(dirrs):
    logfile = dirr+filename
    
    prop_1 = 'Time'
    prop_2 = 'c_PAO_AO_interaction'
    
    thermo = lfp.thermo_panda(logfile, serial=3, zero_ref='Time')
        
    x, y = thermo[prop_1], thermo[prop_2]
    
    start = 1
    end   = None
    
    x, y = x.iloc[start:end], y.iloc[start:end]

    ax.plot(x,y,label=labels[i])
plt.xlabel(f'{ne.getlab(prop_1)}')
plt.ylabel('Interaction Energy (Kcal/mol)')
plt.legend()