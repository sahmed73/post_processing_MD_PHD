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
plt.rcParams['font.size']=15

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

df1   = pd.DataFrame(data)
mean = df1.mean()
std  = df1.std()

df2  = pd.DataFrame(data2)
mean2 = df2.mean()
std2  = df2.std()
df1   = pd.DataFrame()
df1['Scavenging Rate (radicals/ns)'] = mean
df1['std']  = std
df1['Group'] = 'First 2 ns'
df1['Antioxidant'] = AOs

df2   = pd.DataFrame()
df2['Scavenging Rate (radicals/ns)'] = mean2
df2['std']  = std2
df2['Group'] = 'Last 2 ns'
df2['Antioxidant'] = AOs

# ax.bar(AOs, mean, yerr=std, capsize=5, color=colors,
#        edgecolor='black', alpha=0.7)
# ax.set_xlabel('Antioxidants')
# ax.set_ylabel('Scavenging Rate (radicals/ns)')
# ax.set_title('Scavenging Rate During the First 2 ns\n', fontsize=14)
#%%

df = pd.concat([df1[['Antioxidant', 'Scavenging Rate (radicals/ns)', 'std', 'Group']], 
                df2[['Antioxidant', 'Scavenging Rate (radicals/ns)', 'std', 'Group']]])

sns.set(style="whitegrid")
fig, ax = plt.subplots(dpi=350)

# Create the bar plot
ax = sns.barplot(x='Antioxidant', y='Scavenging Rate (radicals/ns)', hue='Group', 
                 data=df, capsize=0.2, errwidth=1.5, ci=None, palette='Set2')

# Manually add error bars
for i, (antioxidant, subset) in enumerate(df.groupby('Antioxidant')):
    for j, (_, row) in enumerate(subset.iterrows()):
        x_pos = i + (j - 0.5) * 0.45  # Adjust position for grouped bars
        plt.errorbar(x_pos, row['Scavenging Rate (radicals/ns)'], 
                     yerr=row['std'], fmt='none', capsize=5, 
                     color='black', elinewidth=1.5)

# Labels and title
ax.set_xlabel('Antioxidants')
ax.set_ylabel('Scavenging Rate (radicals/ns)')
ax.set_title('Scavenging Rate During First vs Last 2 ns')

plt.show()
