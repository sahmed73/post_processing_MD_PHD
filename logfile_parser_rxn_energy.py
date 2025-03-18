# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Apr  6 15:26:03 2024
"""
import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne

plt.style.use('default')
plt.rcParams['font.size']=18


dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0001\Reaction\TSD=2.5\Sim-1"
filename = r'\log.lammps'

logfile = dirr+filename

thermo = lfp.thermo_panda(logfile, serial=':')

prop_1 = 'Time'
prop_2 = 'PotEng'
prop_3 = 'v_rxn1'

time, pe, n_rxn = thermo[prop_1]/1000, thermo[prop_2]/1e6, thermo[prop_3]

start = 0
end   = None#1480

time, pe, n_rxn = time.iloc[start:end], pe.iloc[start:end], n_rxn.iloc[start:end]


last_rxn=0
last_pe=0
deltas = []
for i in range(time.size):
    if n_rxn[i]>last_rxn:
        delta=pe[i]-last_pe
        print("delta E =", delta)
        deltas.append(delta)
        
    last_rxn = n_rxn[i]
    last_pe  = pe[i]

fig, ax = plt.subplots(dpi=300)
ax.plot(time, pe,c='tab:blue', label='Potential Energy $\\times 10^{-6}$ (Kcal/mol)')
ax.plot(time, n_rxn,c='tab:red', label='Reaction')


ax.set_xlabel(f'{ne.getlab(prop_1)}')
# ax.set_ylabel(f'{ne.getlab(prop_2)}')
ax.set_xlabel('Time (ns)')
ax.grid(alpha=0.6)
ax.legend(fontsize=12)
plt.show()

plt.plot(deltas, '--', marker='o')







#%%
import numpy as np

# Constants
k_B = 0.001987  # Boltzmann constant in kcal/mol·K
T = 300  # Temperature in Kelvin

# Load your energy data (replace 'pe' with your potential energy column)
energy = thermo['TotEng']  # Replace with the actual file
energy_min = energy.min() 
normalized_energy = energy - energy_min
beta = 1 / (k_B * T)

# Calculate the Boltzmann factors
boltzmann_factors = np.exp(-beta * normalized_energy)

# Partition function
Q = np.sum(boltzmann_factors)

print(f"Partition function (Q): {Q}")
#%% Enthalpy

E = thermo['TotEng'] # total energy
H = E.copy() # NVT ==> the volume is constant

time, pe, n_rxn = thermo[prop_1]/1000, thermo[prop_2], thermo[prop_3]

last_rxn=0
last_pe=0
deltas = []
for i in range(time.size):
    if n_rxn[i]>last_rxn:
        delta=pe[i]-last_pe
        print("delta E =", delta)
        deltas.append(delta)
    last_rxn = n_rxn[i]

delta_E = sum(deltas)