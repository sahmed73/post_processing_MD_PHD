# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Mar 13 06:55:41 2025
"""

import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def first_order_growth(t, k):
    N0 = 50 
    return N0 * (1 - np.exp(-k * t))


# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=18

parent_dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001'

AOs = ['A0001','A0002','A0003','A0004','A0005']
MAP = dict(zip(AOs, ['L1','L2','S1','S2','S3']))
colors = ['tab:blue', 'tab:orange','tab:green','tab:red','tab:purple']
fig, ax = plt.subplots(dpi=300)


for AO in AOs:
    data = pd.DataFrame()
    for i in range(3):
        dirr = parent_dirr+rf'\{AO}\Reaction\TSD=1.95\Sim-{i+1}'
    
        filename = r'\log.lammps'
        logfile = dirr+filename
        
        thermo = lfp.thermo_panda(logfile, serial=1)
        
        prop_1 = 'Time'
        prop_2 = 'v_rxn1'
        
        x, y = thermo[prop_1]/1000, thermo[prop_2]
        if AO in ['A0002','A0005']:
            y+=thermo['v_rxn2']
        
        data[f'Sim-{i+1}'] = y
    y_mean = data.mean(axis=1)
    y_std = data.std(axis=1)
    ax.scatter(x[::3500], y_mean[::3500],
                            marker='.', label=AO)
    
    popt, _ = curve_fit(first_order_growth, x,y)
    y_fit = first_order_growth(x, *popt)
    ax.plot(x,y_fit)
    
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Number of scavenged radical')
ax.set_ylim(-1,21)
ax.legend(fontsize=12, loc='upper left')