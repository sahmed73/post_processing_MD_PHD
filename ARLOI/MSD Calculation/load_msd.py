# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Dec  9 12:33:43 2024
"""
import numpy as np
import matplotlib.pyplot as plt


colors=['tab:blue', 'tab:green']
# Plot the data
fig, ax = plt.subplots(dpi=350)

for i, AO in enumerate(['A0001','A0003']):
    
    data = np.load(f"msd_{AO}.npy")
    lagtimes, msd = data[:, 0], data[:, 1]
    
    # Plot the actual MSD
    ax.plot(lagtimes, msd, label=AO, c=colors[i])

# Add labels and legend
ax.set_xlabel('Time (ns)')
ax.set_ylabel('MSD ($\mathrm{Å}^2$)')
ax.legend()

# Show the plot
plt.show()