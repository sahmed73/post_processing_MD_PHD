# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Sep 27 12:03:23 2024
"""

import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np

def add_dim(text, xy, xytext, ax):
    ax.annotate('', xy=xy, xytext=xytext, 
                arrowprops=dict(arrowstyle='<->', color='black', alpha=0.6, lw=0.8, shrinkA=0, shrinkB=0),
                xycoords='data', textcoords='data')
    ax.text((xy[0]+xytext[0])/2, xy[1]*1.01, text, fontsize=12, ha='center', va='bottom', transform=ax.transData)
    
# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=15

fig, ax = plt.subplots(dpi=350)
AOs = ['A0001']#, 'A0002', 'A0003']
    
for AO in AOs:
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Equilibration\longer_each_steps'.format(AO)
    filename = r'\log.lammps'
    
    # logfile = ucm.local_copy(dirr+filename)
    logfile = dirr+filename
    
    thermo = lfp.thermo_panda(logfile, serial=[1,2,3,4], timestep=0.25,
                              zero_ref='energy')
    
    property = 'Density'
    
    start = 30
    end   = None
    
    x, y = thermo.Time.iloc[start:end], thermo[property].iloc[start:end]
    
    ax.plot(x,y, label=AO, alpha=0.8)
    
ax.set_xlabel(f"{ne.getlab('Time')}")
ax.set_ylabel(f"{ne.getlab(property)}")

y_text=1.09
add_dim('NPT$_{1}$', xy=(0, y_text), xytext=(150, y_text), ax=ax)
add_dim('NPT$_{2}$', xy=(150, y_text), xytext=(200, y_text), ax=ax)
add_dim('NPT$_{3}$', xy=(200, y_text), xytext=(350, y_text), ax=ax)
# add_dim('NVT', xy=(120, y_text), xytext=(220, y_text), ax=ax)
ax.set_ylim(0.2-0.1,1.19)

plt.legend(loc='lower right')
# mplt.saveplot(folder='thermoplot',name='thermo')