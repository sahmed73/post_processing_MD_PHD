# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Dec  9 12:16:06 2024
"""
 
import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np
 
# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=16
 
 
fig, ax = plt.subplots(dpi=300)
tsds={'1.29':1, '2.00':2, '1.40':3, '1.50':4, '1.60':5, '1.70':6, '1.80':7, '1.95':8}
sorted_tsds=sorted(tsds,reverse=True)
 
for tsd in sorted_tsds:
    dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0001\Reaction\Sim-{tsds[tsd]}_TSD={tsd}"
    filename = r'\log.lammps'
    logfile = dirr+filename
    thermo = lfp.thermo_panda(logfile, serial='all', zero_ref='energy+time')
    prop_1 = 'Time'
    prop_2 = 'v_rxn1'
    x, y = thermo[prop_1]/1000, thermo[prop_2]
    start = 0
    end   = 100000
    x, y = x.iloc[start:end], y.iloc[start:end]
    
    ax.plot(x,y, label=tsd)
    
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Number of reactions')
ax.legend(ncol=2, fontsize=13, loc='upper left',title='R$_{max}$')