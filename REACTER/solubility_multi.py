# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Mar  2 20:46:41 2025
"""

import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D  # For custom legend handles

plt.style.use('default')
plt.rcParams['font.size'] = 18
fig, ax = plt.subplots(dpi=350)
markers = ['o', 's', '^']

n_bulk = 100  # Number of molecules in the bulk simulation
molecule = 'Indomethacin'
y_lim = [None, None]

# Store solubility parameters from multiple simulations
solubility_parameters = []

# Define colors for Gas (Blue) and Bulk (Orange)
colors = [(0.121, 0.466, 0.705, 0.4), (1.0, 0.498, 0.055, 0.4)]  # Tab:Blue and Tab:Orange

# Loop over multiple simulations
for i in range(3):
    print('=' * 50)
    sim = f'Sim-{i+1}'
    print(f"Processing {sim}")

    E_values = []

    for j, mode in enumerate(['Gas', 'Bulk']):
        print('-' * 50)
        print(f"{mode} - {sim}")

        dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\Solubility\{molecule}\PCFF-Gastiger\Eq\{mode}_phase\{sim}"
        filename = r'\log.lammps'
        
        logfile = dirr + filename
        thermo = lfp.thermo_panda(logfile, serial=-1, zero_ref='Time')
        
        time, pe = thermo['Time'], thermo['PotEng']
        start = -2000
        end   = None
        ival  = 1

        if mode == 'Bulk':
            pe = pe / n_bulk
            V_bulk = thermo['Volume'].iloc[start:end:ival].mean()  # Average last 1000 values
            eq_density = thermo['Density'].iloc[start:end:ival].mean()  # Average density
            
        time, pe = time.iloc[start:end:ival] - time.iloc[start:end:ival].min(), pe.iloc[start:end:ival]
        E_values.append(pe.mean())

        # Scatter plot with colors and markers
        ax.scatter(time, pe, s=30, label=f"{mode} - {sim}",
                   facecolors=colors[j],  # Sets the face (fill) color
                   edgecolors='k',        # Black edge color
                   linewidths=0.5,        # Bold edge lines
                   marker=markers[i])     # Different marker per simulation
        
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('P.E. per molecule (Kcal/mol)')
        # ax.set_ylim(y_lim)

    # Compute solubility parameter for this simulation
    E_gas, E_bulk = E_values  # kcal/mol
    V_m = (V_bulk / n_bulk) * (6.022e23 * 1e-24)  # cm³ per mole

    conversion_factor = 64.6838  # (MPa^0.5) / (kcal/mol)^0.5 / cm^(3/2)
    delta_E_vap = (E_gas - E_bulk) 
    delta = np.sqrt(delta_E_vap / V_m) * conversion_factor

    # Store solubility parameter from this simulation
    solubility_parameters.append(delta)

    print('-' * 50)
    print(f"Simulation: {sim}")
    print(f"E_gas   = {E_gas:.4f} kcal/mol")
    print(f"E_bulk  = {E_bulk:.4f} kcal/mol")
    print(f"E_vap   = {E_gas - E_bulk:.4f} kcal/mol")

    print(f"V_bulk  = {V_bulk:.4f} Å³")
    print(f"V_m     = {V_m:.4f} cm³/mol")
    print(f"Density = {eq_density:.2f}")

    print(f"Solubility Parameter (δ): {delta:.2f} MPa^0.5")

# Compute the average solubility parameter
avg_delta = np.mean(solubility_parameters)
std_delta = np.std(solubility_parameters)

print('=' * 50)
print(f"Average Solubility Parameter (δ): {avg_delta:.2f} ± {std_delta:.2f} MPa^0.5")

# Create separate legends for colors (Gas vs. Bulk) and markers (Sim-1, Sim-2, Sim-3)
color_legend_handles = [Line2D([0], [0], marker='.', color='w', markerfacecolor=colors[0], markersize=10, label="Gas"),
                        Line2D([0], [0], marker='.', color='w', markerfacecolor=colors[1], markersize=10, label="Bulk")]

marker_legend_handles = [
    Line2D([0], [0], marker=markers[0], color='k', markerfacecolor='none', markersize=10, linestyle='None', label="Sim-1"),
    Line2D([0], [0], marker=markers[1], color='k', markerfacecolor='none', markersize=10, linestyle='None', label="Sim-2"),
    Line2D([0], [0], marker=markers[2], color='k', markerfacecolor='none', markersize=10, linestyle='None', label="Sim-3")
]

# Display two separate legends
legend1 = ax.legend(handles=color_legend_handles, title="Phase", loc="upper right", fontsize=12, title_fontsize=12)
legend2 = ax.legend(handles=marker_legend_handles, title="Simulation", loc="upper left", fontsize=12, title_fontsize=12)

# Add the first legend back to the plot
ax.add_artist(legend1)

plt.show()
