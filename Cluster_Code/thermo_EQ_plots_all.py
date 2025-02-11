# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Sep 15 19:53:23 2024
"""

import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np
import magnolia.access_ucm_cluster as ucm
import os
import sys

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=20

fig, ax = plt.subplots(dpi=350, figsize=[10,6])

ssh = ucm.connect()
remote_dirr = '/mnt/borgstore/amartini/sahmed73/ARLOI-V1/5_ARLOI_SpeedEQ/20_PAO-OH_15_A0005/Equilibrate'

stdin, stdout, stderr = ssh.exec_command(f'find {remote_dirr} -maxdepth 1 -type d')
directories = stdout.read().decode().splitlines()[1:-1]
ssh.close()


for folder in directories:
    filename = '/log.lammps'
    
    sim_name = folder[folder.rfind('/')+1:]
    print(sim_name)
    
    logfile = ucm.local_copy(folder+filename)
    
    thermo = lfp.thermo_panda(logfile, serial=2, timestep=0.25)
    
    prop_1 = 'Time'
    prop_2 = 'Density'
    x, y = thermo[prop_1], thermo[prop_2]
    
    prop_list = ['PotEng','TotEng','KinEng', 'Time']
    if prop_1 in prop_list:  x = x  - x.min()
    if prop_2 in prop_list:  y = y  - y.min()
    
    target_value = 10000
    closest_index = (np.abs(x - target_value)).idxmin()
    
    start = 10
    end   = closest_index
    
    x, y = x.iloc[start:end], y.iloc[start:end]
    
    ax.plot(x,y, label=sim_name)
    
plt.xlabel(f'{ne.getlab(prop_1)}')
plt.ylabel(f'{ne.getlab(prop_2)}')
plt.legend()
mplt.saveplot(folder='thermoplot',name='thermo')