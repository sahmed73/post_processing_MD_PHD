# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Apr  4 13:08:22 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use('classic')
plt_width = 7
aspect_ratio = 1.333333
plt.figure(figsize=[plt_width,plt_width/aspect_ratio])
plt.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.labelpad'] = 15
plt.rc('axes', grid=True)
plt.rc('font', size=10) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=20) # Axes label size
plt.rc('xtick', labelsize=20) # X-axis tick label size
plt.rc('ytick', labelsize=20) # Y-axis tick label size
plt.rc('legend', fontsize=20) # Legend fontsize

# Read data from file into a NumPy array
data = np.loadtxt("temp.txt")

# Plot the data
plt.plot(np.linspace(0,1125,len(data)),data)
plt.xlabel('Time (ps)')
plt.ylabel('HAT Count')
plt.minorticks_on()
plt.savefig('plot.png',dpi=300,transparent=True, bbox_inches='tight')