# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Apr  6 15:26:03 2024
"""
import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('default')
plt.rcParams['font.size']=14
fig, ax = plt.subplots(dpi=300)
colors = [
    "tab:blue", "tab:red", "tab:green", "tab:orange", 
    "tab:purple", "tab:brown", "tab:pink", "tab:gray", 
    "tab:olive", "tab:cyan"
]



temps = np.arange(300,1001, step=100)
sims = ['Sim-6_TSD=1.90', 'Sim-7_T=400K', 'Sim-8_T=500K', 'Sim-9_T=600K',
        'Sim-10_T=700K', 'Sim-11_T=800K', 'Sim-12_T=900K', 'Sim-13_T=1000K']

rates = []

for sim, temp, color in zip(sims, temps, colors):
    dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0003\Reaction\{sim}"
    filename = r'\log.lammps'
    
    logfile = dirr+filename
    
    thermo = lfp.thermo_panda(logfile, serial='all')
    
    prop_1 = 'Time'
    prop_2 = 'v_rxn1'
    x, y = thermo[prop_1]/1000, thermo[prop_2]
    
    m, c = np.polyfit(x, y, 1)
    rates.append(m)
    ax.plot(x,y, label=f'{temp}K', color=color)
    # ax.plot(x, m*x+c, '--', alpha=0.5)

rates=np.array(rates)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Number of reactions')
ax.legend(title='Temperature', fontsize=12, title_fontsize=12, ncol=2)

#%%
fig, ax = plt.subplots(dpi=300)
# ax.plot(temps,rates, linestyle='--', marker='o', color='maroon')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Reaction rates ($ns^{-1}$)') 
ax.set_xticks(temps)


slope, lnA = np.polyfit(1/temps, np.log(rates), 1)
R = 0.001987  # kcal/mol·K
Ea = -slope*R

print(Ea)

#%%
# Constants in kcal/mol·K
k_B = 1.987e-3   # Boltzmann constant (kcal/mol·K)
h = 3.8088e-14   # Planck's constant (kcal·s/mol)

# Compute Q^‡ / Q_reactants
Q_ratio = (rates*1e9  * h) / (k_B * temps) * np.exp(Ea / (k_B * temps))

print(Q_ratio)

# fig, ax = plt.subplots(dpi=300)
ax.plot(temps,Q_ratio*1e4, linestyle='--', marker='o', color='green')
ax.set_ylabel('Q-Ratio') 






