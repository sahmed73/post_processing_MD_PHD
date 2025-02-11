# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  8 01:56:23 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.metrics import r2_score

plt.rcParams['font.size'] = 13
ramp_rate = 1  # Temperature change per unit time
initial_temp = 300  # Starting temperature
ref_time = 100

color_array = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

fig, ax = plt.subplots(dpi=350)

AOs = ['A0001', 'A0003']#, 'A0004', 'A0005']
for i, AO in enumerate(AOs):
    color = color_array[i % len(color_array)]  
    
    dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\13_ARLOI_Fixed-300K\{AO}\Production\300-500K_TRR=1Kpps"
    df = pd.read_csv(dirr + "\\number_of_SR_vs_time.csv", index_col=0)
    
    # Primary plot on the original axis
    mplt.mean_range_plot(df, ax=ax, color=color, label=f'{AO}', shade=False)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Number of Reaction")
    
    # Create a secondary x-axis for temperature
    # secax = ax.secondary_xaxis('top', functions=(lambda x: x * ramp_rate + initial_temp, 
    #                                              lambda x: (x - initial_temp) / ramp_rate))
    # secax.set_xlabel("Temperature (K)")
    # ax.grid(alpha=0.2)
    
    # # Enable minor ticks and set tick parameters for all axes
    # ax.minorticks_on()
    # ax.tick_params(axis='both', which='major', length=6, width=1.5, direction='in', top=True, right=True)
    # ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in', top=True, right=True)
    # secax.tick_params(axis='x', which='major', length=6, width=1.5, direction='in', top=True)
    # secax.tick_params(axis='x', which='minor', length=3, width=1, direction='in', top=True)
    
    # Fitting
    # time = df.index
    # number_of_sr = df.mean(axis=1)
    # coefficients = np.polyfit(time, number_of_sr, 2)
    # quadratic_fit = np.poly1d(coefficients)
    
    # x_fit = np.linspace(time.min(), time.max(), 100)
    # y_fit = quadratic_fit(x_fit)
    
    # # Plotting actual data and quadratic fit with the same color for each AO
    # ax.plot(x_fit, y_fit, color=color, linestyle='--', alpha=0.8)
    # ax.set_ylim(-1, 16)
    # ax.set_xlim(-10, 210)
    
    # x_value = ref_time
    # y_value = quadratic_fit(x_value)
    # ax.plot([x_value,x_value], [0-1,y_value], '--', color='k', alpha=0.5)
    # ax.plot([-10, x_value], [y_value, y_value], '--', color=color, alpha=0.8)
    # ax.plot(x_value, y_value, 'o', color=color, alpha=0.5)
    # print(AO, y_value)
    
    # # R-squared
    # y_fit_values = quadratic_fit(time)
    # r_squared = r2_score(number_of_sr, y_fit_values)

ax.legend(title='Antioxidants')
# ax.grid(alpha=0.3)
plt.show()
