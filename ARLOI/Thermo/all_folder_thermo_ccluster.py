# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Sep 18 00:22:44 2024
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

# Execute the 'ls' command to list directories
ssh = ucm.connect()

# Read the command output
remote_dirr = '/mnt/borgstore/amartini/sahmed73/ARLOI-V1/5_ARLOI_SpeedEQ/20_PAO-OH_15_A0005/Equilibrate'
stdin, stdout, stderr = ssh.exec_command(f'ls -d {remote_dirr}/*/')
directories = stdout.read().decode().splitlines()


fig, ax = plt.subplots(dpi=350, figsize=[10,6])

for dirr in directories:
    label = dirr[dirr.rfind('Equilibrate/')+len('Equilibrate/'):-1]
    print(label)
    if label=='T-1000_P-1':
        
        continue
    
    filename = '/log.lammps'
    logfile = ucm.local_copy(dirr+filename)
    
    thermo = lfp.thermo_panda(logfile, serial=2, timestep=0.25)
    
    prop_1 = 'Time'
    prop_2 = 'PotEng'
    x, y = thermo[prop_1], thermo[prop_2]
    
    prop_list = ['PotEng','TotEng','KinEng', 'Time']
    if prop_1 in prop_list:  x = x  - x.min()
    if prop_2 in prop_list:  y = y  - y.min()
    
    target_value = 10000
    closest_index = (np.abs(x - target_value)).idxmin()
    
    start = 1
    end   = closest_index
    
    x, y = x.iloc[start:end], y.iloc[start:end]
    
    ax.plot(x,y, label=f'{label}')
    
plt.xlabel(f'{ne.getlab(prop_1)}')
plt.ylabel(f'{ne.getlab(prop_2)}')
plt.legend()
mplt.saveplot(folder='thermoplot',name='thermo')