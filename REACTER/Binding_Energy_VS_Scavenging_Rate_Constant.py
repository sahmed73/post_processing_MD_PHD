# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Mar 12 11:14:13 2025
"""
import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=18



AOs = ['A0001','A0002','A0003','A0004','A0005']
MAP = dict(zip(AOs, ['L1','L2','S1','S2','S3']))
colors = ['tab:blue', 'tab:orange','tab:green','tab:red','tab:purple']

fig, ax = plt.subplots(dpi=300)

Avg_Binding_Energy = []

for i, AO in enumerate(AOs):
    
    dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\{AO}\Eq\Interaction_Energy+MSD"

    filename = r'\log.lammps'
    logfile = dirr+filename
    
    thermo = lfp.thermo_panda(logfile, serial=-1) # last serial
    
    prop_1 = 'Time'
    prop_2 = 'c_PAO_AO_interaction'
    
    x, y = thermo[prop_1]/1000, thermo[prop_2]
    
    start = 2
    end   = None
    ival  = 1
    
    x, y = x.iloc[start:end:ival]-x.iloc[start:end:ival].min(), y.iloc[start:end:ival]
    ax.plot(x,-y, label=MAP[AOs[i]])
    
    Avg_Binding_Energy.append(-y.mean())
    # ax.plot(x,[-y.mean()]*x.size,color='k')
    
ax.set_xlabel('Time (ns)')
ax.set_ylabel('$E_{binding}$ (Kcal/mol)')
# ax.legend(fontsize=12)
ax.set_ylim(bottom=1001, top=1900)
legend_handles = [Line2D([0], [0], color=colors[i], lw=3.5, label=MAP[AOs[i]]) for i in range(len(AOs))]
ax.legend(handles=legend_handles, fontsize=16)

#%%

fig, ax2 = plt.subplots(dpi=350)
data = {
    "Sim-1": [0.081444, 0.077283, 0.074504, 0.080203, 0.116402],
    "Sim-2": [0.076825, 0.122914, 0.046360, 0.045120, 0.082771],
    "Sim-3": [0.044220, 0.098843, 0.058524, 0.077844, 0.096224]
}

df = pd.DataFrame(data, index=AOs)
k = df.mean(axis=1)
E_int = np.array(Avg_Binding_Energy)

# Apply glowing effect using multiple scatter layers
for alpha, size in zip([0.1, 0.2, 0.3, 0.4, 0.5], [300, 250, 200, 150, 100]):
    ax2.scatter(k, E_int, s=size, color=colors, alpha=alpha,
                edgecolor='k')

ax2.set_xlabel("Scavenging Rate Constant (ns$^{-1}$)")
ax2.set_ylabel("Average $E_{binding}$ (Kcal/mol)")
ax2.set_xlim(left=0.055, right=0.105)
ax2.set_ylim(bottom=1001, top=1900)
ax2.grid(alpha=0.4, linestyle='--')

for i, txt in enumerate(AOs):
    txt = MAP[txt]
    if txt=='S3': offset = (-15,-23)
    else: offset = (5,5)
    ax2.annotate(txt, 
                 (k[i], E_int[i]),  # Position
                 fontsize=14, 
                 ha='left', va='bottom',  # Adjust alignment
                 xytext=offset,  # Offset text slightly (x, y)
                 textcoords='offset points')

