# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 22 10:46:59 2024
"""

import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=20


dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\TEST\T1\Geometry\bulk-03\REACTION-v3'
filename = r'\reaction.log.lammps'

logfile = dirr+filename

thermo = lfp.thermo_panda(logfile, serial=1, timestep=0.5,
                          zero_ref='energy+time')

prop_1 = 'Time'
prop_2 = 'v_rxn1'
x, y = thermo[prop_1]/1000, thermo[prop_2]

start = 0
end   = None

x, y = x.iloc[start:end], y.iloc[start:end]

fig, ax = plt.subplots(dpi=300)
ax.plot(x,y,c='tab:red')
# plt.xlabel(f'{ne.getlab(prop_1)}')
# plt.ylabel(f'{ne.getlab(prop_2)}')
plt.xlabel('Time (ns)')
plt.ylabel('Number of reactions')

def weibull_func(t, L, t0, k):
    return L * (1 - np.exp(-(t / t0) ** k))



# Perform curve fitting
params, covariance = curve_fit(weibull_func, x, y, p0=(1, 0.1, 1),
                               maxfev=5000)

# Generate fitted y values
y_fit = weibull_func(x, *params)

# Plot the data and the fitted curve
fig, ax = plt.subplots(dpi=300)
ax.plot(x, y, label='Data', c='tab:blue')
ax.plot(x, y_fit, label='Fit', c='tab:red', linestyle='--')

# Labels and legend
plt.xlabel('Time (ns)')
plt.ylabel('Number of reactions')
plt.legend()

# Save the plot
plt.savefig('thermo_with_fit.png', dpi=300)

# Print fit parameters
print('Fit parameters (a, b, c):', params)
