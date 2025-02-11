# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Sep 15 20:37:18 2024
"""
import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np
import magnolia.access_ucm_cluster as ucm

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=20

AOs = ['A0006', 'A0007', 'A0008', 'A0009', 'A0010']
colors=[]
fig, ax = plt.subplots(dpi=350, figsize=[10,6])

for i, AO in enumerate(AOs):
    dirr = '/mnt/borgstore/amartini/sahmed73/ARLOI-V2/1_ARLOI_24-09-12--00-35-24/20_PAO-OH_15_{}/Equilibrate/Sim-2'.format(AO)
    filename = '/log.lammps'
    logfile = ucm.local_copy(dirr+filename)
    
    thermo = lfp.thermo_panda(logfile, serial=3, timestep=0.25)
    
    prop_1 = 'Time'
    prop_2 = 'Press'
    x, y = thermo[prop_1], thermo[prop_2]
    
    prop_list = ['PotEng','TotEng','KinEng', 'Time']
    if prop_1 in prop_list:  x = x  - x.min()
    if prop_2 in prop_list:  y = y  - y.min()
    
    target_value = 10000
    closest_index = (np.abs(x - target_value)).idxmin()
    
    start = 10
    end   = closest_index
    
    x, y = x.iloc[start:end], y.iloc[start:end]
    
    ax.plot(x,y, label=f'{AO}')
    plt.xlabel(f'{ne.getlab(prop_1)}')
    plt.ylabel(f'{ne.getlab(prop_2)}')
    plt.legend()
    mplt.saveplot(folder='thermoplot',name='thermo')