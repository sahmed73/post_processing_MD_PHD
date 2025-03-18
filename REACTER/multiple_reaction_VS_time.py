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
import pandas as pd

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=23

parent_dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001'

AOs = ['A0001','A0002','A0003','A0004','A0005']
MAP = dict(zip(AOs, ['L1','L2','S1','S2','S3']))
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
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
            # y = y/2
        
        data[f'Sim-{i+1}'] = y
    data.index=x
    # mplt.mean_range_plot(data, ax=ax, shade=True,
    #                      label=AO, alpha=0.2)
    
    ### errorbar
    x = data.index
    y_mean = data.mean(axis=1)
    y_std = data.std(axis=1)
    _, caps, bars = ax.errorbar(x[::3500], y_mean[::3500], y_std[::3500],
                            marker='.', capsize=3, label=MAP[AO], 
                            elinewidth=1)  # Set color explicitly if needed

    # Adjust transparency only for error bars (bars and caps)
    for bar in bars:
        bar.set_alpha(0.0)  # Make error bars semi-transparent
    for cap in caps:
        cap.set_alpha(0.0)  # Make caps semi-transparent
    
    ## scatter plot
    # x = data.index
    # y_mean = data.mean(axis=1)
    # ax.plot(x[::1], y_mean[::1],
    #         marker='.', markevery=3500, label=AO)
    
ax.set_xlabel('Time (ns)')
ax.set_ylabel('# Radicals scavenged')
ax.set_ylim(-1,21)
ax.set_xticks(range(6))
ax.legend(fontsize=18, loc='upper left')