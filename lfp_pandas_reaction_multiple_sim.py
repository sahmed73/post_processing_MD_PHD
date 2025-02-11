# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Dec  7 22:11:41 2024
"""

import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=16


fig, ax = plt.subplots(dpi=300)
tsds=['1.80', '1.95']
AOs={'A0001':7, 'A0003':6}

for AO, tsd in zip(AOs,tsds):
    print(AO,tsd)
    dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\{AO}\Reaction\Sim-{AOs[AO]}_TSD={tsd}"
    filename = r'\log.lammps'
    logfile = dirr+filename
    thermo = lfp.thermo_panda(logfile, serial='all', zero_ref='energy+time')
    prop_1 = 'Time'
    prop_2 = 'v_rxn1'
    x, y = thermo[prop_1]/1000, thermo[prop_2]
    start = 0
    end   = 100000
    x, y = x.iloc[start:end], y.iloc[start:end]
    
    ax.plot(x,y, label=AO)
    
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Number of reactions')
ax.legend(ncol=1, fontsize=13, loc='upper left',title='Antioxidants')