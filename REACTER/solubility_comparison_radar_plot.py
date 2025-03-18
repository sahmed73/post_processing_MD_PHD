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
reported_values = [34.8, 23.9, 13.50]  # Reported MD Sim Values
# experimental_values = [30.0, 22.0, 10.0]  # Assumed experimental MPa^0.5 values

# Number of variables
num_vars = len(molecules)

# Compute angles for radar chart
angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

# Close the radar chart for all data sets
simulated_values += [simulated_values[0]]
reported_values += [reported_values[0]]
# experimental_values += [experimental_values[0]]
angles += [angles[0]]

# Create radar plot
fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

# Plot simulated data
ax.fill(angles, simulated_values, color='blue', alpha=0.3, label="Simulated")
ax.plot(angles, simulated_values, color='blue', linewidth=2)

# Plot reported data
ax.fill(angles, reported_values, color='red', alpha=0.3, label="Reported")
ax.plot(angles, reported_values, color='red', linewidth=2)

# Plot experimental data
# ax.fill(angles, experimental_values, color='green', alpha=0.3, label="Experimental")
# ax.plot(angles, experimental_values, color='green', linewidth=2)

# Labels
ax.set_xticks(angles[:-1])
ax.set_xticklabels(molecules, fontsize=14, fontweight='bold', ha='center', va='top')

# Adjust labels to avoid overlap
for label, angle in zip(ax.get_xticklabels(), angles[:-1]):
    label.set_rotation(np.degrees(angle) + 10)  # Slight rotation for better readability
    label.set_verticalalignment("bottom" if np.pi/2 <= angle <= 3*np.pi/2 else "top")

# Title and legend
plt.title("Solubility Parameter Comparison (MPa$^{0.5}$)\n\n", fontsize=14)
plt.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))

# Show plot
plt.show()