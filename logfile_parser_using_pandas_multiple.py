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
plt.rcParams['font.size']=20


dirrs = [r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0001\Eq\Sim-2",
        r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0003\Eq\Sim-1"]

AOs = ['A0001', 'A0003']

fig, ax = plt.subplots(dpi=300)
for dirr, AO in zip(dirrs,AOs):        
    filename = r'\log.lammps'
    
    logfile = dirr+filename
    
    thermo = lfp.thermo_panda(logfile, serial='all', zero_ref='energy+time')
    
    prop_1 = 'Time'
    prop_2 = 'Density'
    x, y = thermo[prop_1]/1000, thermo[prop_2]
    
    start = 0
    end   = None

    x, y = x.iloc[start:end], y.iloc[start:end]
    ax.plot(x,y,label=AO)
    
plt.xlabel(f'{ne.getlab(prop_1)}')
plt.ylabel(f'{ne.getlab(prop_2)}')
plt.xlabel('Time (ns)')
plt.legend()
# plt.ylabel('Number of reactions')
mplt.saveplot(folder='thermoplot',name='thermo')
# print('Final Density:',thermo['Density'].iloc[-1])