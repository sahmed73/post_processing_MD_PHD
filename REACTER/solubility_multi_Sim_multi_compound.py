# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Mar 10 10:06:31 2025
"""

import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D  

def count_molecules_and_total_mass(lammps_file):
    """
    Count the number of unique molecules and calculate total mass from a LAMMPS data file.

    Parameters:
        lammps_file (str): Path to the LAMMPS data file.

    Returns:
        tuple: (Number of unique molecules, Total mass in atomic mass units (amu))
    """
    with open(lammps_file, 'r') as file:
        lines = file.readlines()

    # Extract Atomic Masses
    mass_section_start = None
    atom_masses = {}
    
    for i, line in enumerate(lines):
        if line.strip().startswith("Masses"):
            mass_section_start = i + 2  
            break

    if mass_section_start is not None:
        for line in lines[mass_section_start:]:
            if line.strip() == "":
                break  
            parts = line.split()
            atom_type = int(parts[0])  
            mass = float(parts[1])     
            atom_masses[atom_type] = mass

    if not atom_masses:
        raise ValueError("No 'Masses' section found in the LAMMPS data file.")

    # Find the "Atoms" Section
    atom_section_start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atom_section_start = i + 2  
            break

    if atom_section_start is None:
        raise ValueError("No 'Atoms' section found in the LAMMPS data file.")

    # Read Atom Data and Compute Total Mass
    atom_data = []
    total_mass = 0

    for line in lines[atom_section_start:]:
        if line.strip() == "":
            break  
        atom_info = line.split()
        molecule_id = int(atom_info[1])
        atom_type = int(atom_info[2])  
        total_mass += atom_masses[atom_type]  
        atom_data.append(molecule_id)

    # Count Unique Molecules
    unique_molecules = len(set(atom_data))

    return unique_molecules, total_mass


plt.style.use('default')
plt.rcParams['font.size'] = 12

markers = ['o', 's', '^']

# Store solubility parameters from multiple simulations
solubility_parameters = []
avg_deltas = []
std_deltas = []
E_gas_values = []
E_bulk_values = []
V_m_values = []
avg_density = []
std_density = []

AOs = ['PAO', 'A0001', 'A0002', 'A0003', 'A0004', 'A0005']
# AOs = ['Chloroethane','Glucose','Indomethacin']
colors = ['grey', 'tab:blue', 'tab:orange','tab:green','tab:red','tab:purple']

for molecule in AOs:
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
                min_data = dirr + r'\min.data'
                n_bulk, total_mass = count_molecules_and_total_mass(min_data)  
                
                density = thermo['Density'].iloc[start:end:ival].mean()
                V_occ = total_mass / density  
                densities.append(density)
                
                print()
                print(f"[DEBUG] {molecule} - {sim}")
                print(f"  - n_bulk = {n_bulk}")
                print(f"  - Total Mass = {total_mass}")
                print(f"  - Density = {density}")
                print(f"  - V_occ = {V_occ}")                
                
                pe = pe / n_bulk
                V_bulk = thermo['Volume'].iloc[start:end:ival].mean()
                print(f"  - V_bulk = {V_bulk}")
            
            time, pe = time.iloc[start:end:ival] - time.iloc[start:end:ival].min(), pe.iloc[start:end:ival]
            E_values.append(pe.mean())

        # Energy values
        E_gas, E_bulk = E_values  
        V_m = (V_bulk / n_bulk) * (6.022e23 * 1e-24)  

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
AO_Sol = np.array(avg_deltas[1:])-avg_deltas[0]
fig, ax2 = plt.subplots(dpi=350)
ax2.bar(AOs[1:], AO_Sol, capsize=5, color=colors[1:], edgecolor='black')
ax2.set_ylabel('HSP (MPa$^{1/2}$)')
ax2.set_xlabel('Antioxidant')

plt.show()
#%%


