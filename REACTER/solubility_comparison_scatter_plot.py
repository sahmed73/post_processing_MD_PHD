# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Mar  3 04:39:47 2025
"""

import numpy as np
import matplotlib.pyplot as plt

# Data for solubility (MPa^0.5)
molecules = ["Glucose", "Indomethacin", "Chloroethane"]
simulated_values = [29.68, 25.79, 10.54]  # My MD Sim Values
reported_values = [34.8, 23.9, 13.50]

# Error estimates (assumed for demonstration)
simulated_errors = [0.45, 0.80, 1.05]
reported_errors = [0.2, 0.3, 0.41]

# Define colors for each molecule
colors = ['tab:blue', 'tab:orange', 'tab:green']

# Create scatter plot
fig, ax = plt.subplots(figsize=(6, 6))

# Error bars only
for i in range(len(molecules)):
    ax.errorbar(simulated_values[i], reported_values[i], xerr=simulated_errors[i], yerr=reported_errors[i],
                fmt='none', capsize=5, color='black', alpha=0.6)
    ax.scatter(simulated_values[i], reported_values[i], facecolors=colors[i], edgecolors='k', 
               s=80, linewidths=1.5, label=molecules[i])

# Plot y = x line
min_val = min(simulated_values + reported_values) * 0.9
max_val = max(simulated_values + reported_values) * 1.1
ax.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="black")

# Labels and formatting
ax.set_xlabel("Simulated Solubility Parameter (MPa$^{0.5}$)", fontsize=12)
ax.set_ylabel("Reported Solubility Parameter (MPa$^{0.5}$)", fontsize=12)
ax.set_title("Simulated vs Reported Solubility Parameters", fontsize=14)
ax.legend(loc='lower right')
ax.grid(True, linestyle='--', alpha=0.5)

# Show plot
plt.show()
