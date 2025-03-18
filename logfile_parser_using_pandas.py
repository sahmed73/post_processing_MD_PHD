# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Apr  6 15:26:03 2024
"""
import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np
import pandas as pd

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=18


dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0001\Reaction\TSD=2.5\Sim-1"
filename = r'\log.lammps'

logfile = dirr+filename

thermo = lfp.thermo_panda(logfile, serial=':', zero_ref='Time')

prop_1 = 'Time'
prop_2 = 'v_rxn1' #'c_PAO_AO_interaction'
x, y = thermo[prop_1], thermo[prop_2]

start = 0 # -2000
end   = -1
ival  = 10 

x, y = x.iloc[start:end:ival]/1000, y.iloc[start:end:ival]

fig, ax = plt.subplots(dpi=350)
ax.plot(x,y,c='tab:red')
# ax.scatter(x,y,c='tab:red', edgecolor='k', linewidth=0.5, alpha=0.7, s=10)
# ax.plot(x, y.rolling(window=50).mean(), color='b')

# ax.plot([1000,1300], [y.mean(), y.mean()])
# plt.xlabel(f'{ne.getlab(prop_1)}')
plt.ylabel(f'{ne.getlab(prop_2)}')

# print(f"Stable Density: {thermo['Density'].iloc[-100:].mean()}")
plt.xlabel('Time (ns)')
# plt.ylabel('Number of reactions')
