# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Feb 23 13:25:15 2025
"""

import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

plt.style.use('default')

parent_dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001'

AOs = ['A0001','A0002','A0003','A0004','A0005']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

data  = {}
data2 = {}
for AO in AOs:
    data[AO] = {}
    data2[AO] = {}
    for i in range(3):
        dirr = parent_dirr+rf'\{AO}\Reaction\TSD=1.95\Sim-{i+1}'
        filename = r'\log.lammps'
        logfile = dirr+filename
        thermo = lfp.thermo_panda(logfile, serial=1)

        prop_1 = 'Time'
        prop_2 = 'v_rxn1'
    

        x, y = thermo[prop_1]/1000, thermo[prop_2]
        if AO in ['A0002','A0005']:
            y += thermo['v_rxn2']
        
        # 1st 2 ns
        slope, _ = np.polyfit(x[:40000], y[:40000], 1)
        data[AO][f'Sim-{i+1}'] = slope
        
        # last 2 ns
        slope, _ = np.polyfit(x[-40000:], y[-40000:], 1)
        data2[AO][f'Sim-{i+1}'] = slope

df1   = pd.DataFrame(data).T
df1['Group'] = 'First 2ns'
df2  = pd.DataFrame(data2).T
df2['Group'] = 'Last 2ns'

df = pd.concat([df1, df2]).reset_index()
df = df.melt(id_vars=['index', 'Group'], var_name='Simulation', value_name='Scavenging Rate')
df.rename(columns={'index': 'Antioxidant'}, inplace=True)

sns.set(style="whitegrid")
fig, ax = plt.subplots(dpi=350)
ax = sns.barplot(x='Antioxidant', y='Scavenging Rate', hue='Group', 
                 data=df, capsize=0.2, err_kws={'linewidth': 1.5}, errorbar=('sd'), palette='Set2')

ax.set_xlabel('Antioxidants', fontsize=15)
ax.set_ylabel('Scavenging Rate (radicals/ns)', fontsize=15)
ax.tick_params(axis='both', labelsize=15)
ax.legend(title=None, fontsize=15, title_fontsize=16)

ax.set_ylim(0,5.99)
