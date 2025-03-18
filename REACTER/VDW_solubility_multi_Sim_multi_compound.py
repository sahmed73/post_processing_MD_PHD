# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Mar 10 10:06:31 2025
"""

import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

plt.style.use('default')
plt.rcParams['font.size'] = 18

# Store solubility parameters from multiple simulations
solubility_parameters = []
avg_deltas = []
std_deltas = []
E_gas_values = []
E_bulk_values = []
V_m_values = []
avg_density = []
std_density = []

n_bulk = 50
AOs = ['PAO', 'A0001', 'A0002', 'A0003', 'A0004', 'A0005']
SMILES = ['CCCCCCCCCCC(CCCCCCCC)CC(C)CCCCCCCC',
          'CCOC(=O)CCC1=CC(=C(O)C(=C1)C(C)(C)C)C(C)(C)C',
          'CC(C)(C)C1=CC(\C=C/CO)=CC(=C1O)C1=CC(\C=C/CO)=CC(=C1O)C(C)(C)C',
          'Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1',
          'CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)C=CC(=O)',
          'CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)CC2=CC(=C(C(=C2)C(C)(C)C)O)C(C)(C)C']
MAP = dict(zip(AOs[1:], ['L1','L2','S1','S2','S3']))
colors = ['grey', 'tab:blue', 'tab:orange','tab:green','tab:red','tab:purple']

for molecule, smile in zip(AOs,SMILES):
    print('\n\n')
    print("*"*35, molecule ,"*"*35)
    solubility_parameters = []
    densities = []
    for i in range(3):
        sim = f'Sim-{i+1}'
        E_values = []
        
        for j, mode in enumerate(['Gas', 'Bulk']):
            dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\Solubility\Solubility_Automation\Solubility001\{molecule}\Solubility_Simulations\{mode}\{sim}"
            filename = r'\\log.lammps'
            
            logfile = dirr + filename
            thermo = lfp.thermo_panda(logfile, serial=-1, zero_ref='Time')
            
            time, pe = thermo['Time'], thermo['PotEng']
            start = -2000
            end   = None
            ival  = 1

            if mode == 'Bulk':
                mol = Chem.MolFromSmiles(smile)
                print(smile)
                V_vdw = Descriptors.MolMR(mol) * 1.3  # Convert to Å³

                density = thermo['Density'].iloc[start:end:ival].mean()
                densities.append(density)
                
                print()
                print(f"[DEBUG] {molecule} - {sim}")
                print(f"  - n_bulk = {n_bulk}")
                print(f"  - Density = {density}")
                print(f"  - V_vdw = {V_vdw}")                
                
                pe = pe / n_bulk
                V_bulk = thermo['Volume'].iloc[start:end:ival].mean()
                print(f"  - V_bulk = {V_bulk}")
            
            time, pe = time.iloc[start:end:ival] - time.iloc[start:end:ival].min(), pe.iloc[start:end:ival]
            E_values.append(pe.mean())

        # Energy values
        E_gas, E_bulk = E_values  
        V_m = (V_vdw * 1e-24) * (6.022e23) 

        # Debugging V_m
        print(f"  - V_m = {V_m} cm³/mol")

        E_gas_values.append(E_gas)
        E_bulk_values.append(E_bulk)
        V_m_values.append(V_m)

        conversion_factor = 64.6838  
        delta_E_vap = (E_gas - E_bulk)
        delta = np.sqrt(delta_E_vap / V_m) * conversion_factor
        
        solubility_parameters.append(delta)

        print(f"  - E_gas = {E_gas} kcal/mol")
        print(f"  - E_bulk = {E_bulk} kcal/mol")
        print(f"  - E_vap = {E_gas - E_bulk} kcal/mol")
        print(f"  - HSP = {delta} MPa^0.5")

    avg_delta = np.mean(solubility_parameters)
    std_delta = np.std(solubility_parameters)
    
    avg_density = np.mean(densities)
    std_density = np.std(densities)
    
    print('='*60)
    print(f'Average delta of {molecule}: {avg_delta:.2f} ± {std_delta:.2f}')
    print('='*60)
    print(f'Average density of {molecule}: {avg_density:.2f} ± {std_density:.2f}')
    print('='*60)
    avg_deltas.append(avg_delta)
    std_deltas.append(std_delta)

# Compute E_vap
E_gas_values = np.array(E_gas_values).reshape(len(AOs), 3)
E_bulk_values = np.array(E_bulk_values).reshape(len(AOs), 3)
V_m_values = np.array(V_m_values).reshape(len(AOs), 3)
E_vap_values = E_gas_values - E_bulk_values

E_vap_means = np.mean(E_vap_values, axis=1)
E_vap_stds = np.std(E_vap_values, axis=1)
V_m_means = np.mean(V_m_values, axis=1)
V_m_stds = np.std(V_m_values, axis=1)

# Plot HSP
fig, ax = plt.subplots(dpi=350)
ax.bar(AOs, avg_deltas, yerr=std_deltas, capsize=5, color=colors, edgecolor='black')
ax.set_ylabel('HSP (MPa$^{1/2}$)')
ax.set_xlabel('Antioxidant')

# Plot E_vap
fig, ax1 = plt.subplots(dpi=350)
ax1.bar(AOs, E_vap_means, yerr=E_vap_stds, capsize=5, color=colors, edgecolor='black')
ax1.set_ylabel('E$_{vap}$ (kcal/mol)')
ax1.set_xlabel('Antioxidant')
#%%
# Plot Solubility
AO_Sol = np.abs(np.array(avg_deltas[1:])-avg_deltas[0])
fig, ax2 = plt.subplots(dpi=500)
ax2.bar([MAP[X] for X in AOs[1:]], AO_Sol, capsize=5, color=colors[1:], edgecolor='black')
ax2.set_ylabel(r'$|\delta-\delta_{PAO}|$ (MPa$^{1/2}$)')
ax2.set_xlabel('Antioxidant')
ax2.set_ylim((0,6.5))

plt.show()
#%%
fig, ax2 = plt.subplots(dpi=350)
data = {
    "Sim-1": [0.081444, 0.077283, 0.074504, 0.080203, 0.116402],
    "Sim-2": [0.076825, 0.122914, 0.046360, 0.045120, 0.082771],
    "Sim-3": [0.044220, 0.098843, 0.058524, 0.077844, 0.096224]
}

df = pd.DataFrame(data, index=AOs[1:])
k = df.mean(axis=1)
E_coh = AO_Sol

# Apply glowing effect using multiple scatter layers
for alpha, size in zip([0.1, 0.2, 0.3, 0.4, 0.5], [300, 250, 200, 150, 100]):
    ax2.scatter(k, E_coh, s=size, color=colors[1:], alpha=alpha,
                edgecolor='k')

ax2.set_xlabel("Scavenging Rate Constant (ns$^{-1}$)")
ax2.set_ylabel(r'$|\delta-\delta_{PAO}|$ (MPa$^{1/2}$)')
ax2.set_xlim(left=0.055, right=0.105)
ax2.set_ylim((0,6.5))
ax2.grid(alpha=0.4, linestyle='--')

for i, txt in enumerate(AOs[1:]):
    txt = MAP[txt]
    if txt in []: offset = (-8,-25)
    else: offset = (5,5)
    ax2.annotate(txt, 
                 (k[i], E_coh[i]),  # Position
                 fontsize=14, 
                 ha='left', va='bottom',  # Adjust alignment
                 xytext=offset,  # Offset text slightly (x, y)
                 textcoords='offset points')

# Perform linear regression
from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(k, E_coh)
r_squared = r_value**2
print(f"R² : {r_squared:.3f}")
# ax2.text(0.05, 0.95, f"R² = {r_squared:.3f}", transform=ax2.transAxes, fontsize=12,
#          verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))



from scipy.stats import pearsonr
corr, _ = pearsonr(k, E_coh)
print(f"Pearson Correlation: {corr:.3f}")

