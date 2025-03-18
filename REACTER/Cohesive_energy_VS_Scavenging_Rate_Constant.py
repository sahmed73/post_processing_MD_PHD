# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Mar 10 10:06:31 2025
"""

import magnolia.log_parser_SA as lfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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


# List of antioxidants
AOs = ['A0001', 'A0002', 'A0003', 'A0004', 'A0005']
MAP = dict(zip(AOs, ['L1','L2','S1','S2','S3']))
colors = ['tab:blue', 'tab:orange','tab:green','tab:red','tab:purple']

# Store cohesive energy values
E_gas_values = []
E_bulk_values = []
cohesive_energy_values = []

for molecule in AOs:
    print("\n" + "*" * 35, molecule, "*" * 35)

    E_gas_list = []  # Store E_gas from different simulations
    E_bulk_list = []  # Store E_bulk from different simulations

    for i in range(3):  # Three simulations per molecule
        sim = f'Sim-{i+1}'
        E_values = []  # Store E_gas and E_bulk for each sim
        
        for mode in ['Gas', 'Bulk']:  # Loop through both Gas and Bulk simulations
            dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\Solubility\Solubility_Automation\Solubility001\{molecule}\Solubility_Simulations\{mode}\{sim}"
            filename = r'\\log.lammps'
            
            logfile = dirr + filename
            thermo = lfp.thermo_panda(logfile, serial=-1, zero_ref='Time')

            pe = thermo['PotEng']
            start = -2000
            E_mean = pe.iloc[start:].mean()  # Average over last 2000 frames

            if mode == 'Bulk':
                min_data = dirr + r'\min.data'
                n_bulk, total_mass = count_molecules_and_total_mass(min_data)  # Get molecule count
                
                E_mean /= n_bulk  # Convert to per-molecule energy
                print(f"  - {sim} | {mode} Energy (per molecule): {E_mean:.2f} kcal/mol")

            else:
                print(f"  - {sim} | {mode} Energy: {E_mean:.2f} kcal/mol")

            E_values.append(E_mean)

        # Store energy values
        E_gas_list.append(E_values[0])  # Gas phase energy
        E_bulk_list.append(E_values[1])  # Bulk phase energy

    # Compute mean and standard deviation per molecule
    E_gas_mean, E_gas_std = np.mean(E_gas_list), np.std(E_gas_list)
    E_bulk_mean, E_bulk_std = np.mean(E_bulk_list), np.std(E_bulk_list)

    # Compute cohesive energy per molecule
    E_cohesive_mean = E_gas_mean - E_bulk_mean
    E_cohesive_std = np.sqrt(E_gas_std**2 + E_bulk_std**2)  # Propagation of uncertainty

    cohesive_energy_values.append((E_cohesive_mean, E_cohesive_std))

    print(f"  >>> {molecule} | Cohesive Energy per molecule: {E_cohesive_mean:.2f} ± {E_cohesive_std:.2f} kcal/mol")

# Convert results to DataFrame
df = pd.DataFrame({
    "Mean Cohesive Energy (per molecule)": [x[0] for x in cohesive_energy_values],
    "Std Dev": [x[1] for x in cohesive_energy_values]
}, index=AOs)

print("\nFinal Cohesive Energy Results (Per Molecule):\n", df)

fig, ax = plt.subplots(dpi=350)
x = [MAP[AO] for AO in AOs]
ax.bar(x, df.iloc[:,0], yerr=df.iloc[:,1],
       color=colors, capsize=5, edgecolor='k', linewidth=1)
ax.set_ylabel('E$_{cohesive}$ (kcal/mol)')
ax.set_xlabel('Antioxidant')

##================
fig, ax2 = plt.subplots(dpi=350)
data = {
    "Sim-1": [0.081444, 0.077283, 0.074504, 0.080203, 0.116402],
    "Sim-2": [0.076825, 0.122914, 0.046360, 0.045120, 0.082771],
    "Sim-3": [0.044220, 0.098843, 0.058524, 0.077844, 0.096224]
}

df = pd.DataFrame(data, index=AOs)
k = df.mean(axis=1)
E_coh = np.array([x for x,y in cohesive_energy_values])

# Apply glowing effect using multiple scatter layers
for alpha, size in zip([0.1, 0.2, 0.3, 0.4, 0.5], [300, 250, 200, 150, 100]):
    ax2.scatter(k, E_coh, s=size, color=colors, alpha=alpha,
                edgecolor='k')

ax2.set_xlabel("Scavenging Rate Constant (ns$^{-1}$)")
ax2.set_ylabel('Average E$_{cohesive}$ (kcal/mol)')
ax2.set_xlim(left=0.055, right=0.105)
ax2.set_ylim(top=31)
ax2.grid(alpha=0.4, linestyle='--')

for i, txt in enumerate(AOs):
    txt = MAP[txt]
    if txt in ['S3','L1']: offset = (-8,-25)
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
