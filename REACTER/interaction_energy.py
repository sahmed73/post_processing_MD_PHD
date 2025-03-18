# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Apr  6 15:26:03 2024
"""
import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=18



AOs = ['A0001','A0002','A0003','A0004','A0005']
fig, ax = plt.subplots(dpi=300)

for i, AO in enumerate(AOs):
    
    dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\{AO}\Eq\Interaction_Energy+MSD"

    filename = r'\log.lammps'
    logfile = dirr+filename
    
    thermo = lfp.thermo_panda(logfile, serial=3)
    
    prop_1 = 'Time'
    prop_2 = 'c_PAO_AO_interaction'
    
    x, y = thermo[prop_1]/1000, thermo[prop_2]
    
    start = 2
    end   = None
    
    x, y = x.iloc[start:end], y.iloc[start:end]
    ax.plot(x,y, label=AOs[i])
    
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Interaction Energy (Kcal/mol)')
ax.legend(fontsize=12)